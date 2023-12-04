// 定义研究区的坐标范围
var roiCoordinates = [
  [122.38155416405702, 39.79603431229974],
  [122.61501363671327, 39.79603431229974],
  [122.61501363671327, 39.99306192199326],
  [122.38155416405702, 39.99306192199326],
  [122.38155416405702, 39.79603431229974]
];
// 创建研究区的多边形
var roi = ee.Geometry.Polygon(roiCoordinates);
Map.addLayer(roi_1, {color: "red"}, "roi_1");
Map.centerObject(roi_1,10);

// remove cloud for Landsat 4, 5 and 7
var rmL457Cloud = function(image) {
    var qa = image.select('pixel_qa');
    // If the cloud bit (5) is set and the cloud confidence (7) is high
    // or the cloud shadow bit is set (3), then it's a bad pixel.
    var cloud = qa.bitwiseAnd(1 << 5)
                    .and(qa.bitwiseAnd(1 << 7))
                    .or(qa.bitwiseAnd(1 << 3));
    // Remove edge pixels that don't occur in all bands
    var mask2 = image.mask().reduce(ee.Reducer.min());
    
    // remove pixels where the blue reflectance is greater than 0.2
    var mask3 = image.select('B1').gt(2000);
    return image.updateMask(cloud.not()).updateMask(mask2).updateMask(mask3.not())
                .copyProperties(image)
                .copyProperties(image, ["system:time_start",'system:time_end','system:footprint']);
  };
  
  // reomove cloud for Landsat-8
  function rmL8Cloud(image) { 
    var cloudShadowBitMask = (1 << 3); 
    var cloudsBitMask = (1 << 5); 
    var qa = image.select('pixel_qa'); 
    var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0) 
                   .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
    var mask2 = image.select('B1').gt(2000);
    return image.updateMask(mask).updateMask(mask2.not())
                .copyProperties(image)
                .copyProperties(image, ["system:time_start",'system:time_end']);
  } 
  var l8_sr = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR").map(rmL8Cloud)
              
              
  var l7_sr = ee.ImageCollection("LANDSAT/LE07/C01/T1_SR").map(rmL457Cloud)
              .filter(ee.Filter.lte('CLOUD_COVER',5))//云量过滤;
              
  var l5_sr = ee.ImageCollection("LANDSAT/LT05/C01/T1_SR").map(rmL457Cloud)
              .filter(ee.Filter.lte('CLOUD_COVER',5))//云量过滤;
              
  
  
  //计算水体方法
  //MNDWI>NDVI
  function calcWater(image) {
    var MNDWI = image.select("MNDWI");
    var NDVI = image.select("NDVI");
    var EVI = image.select("EVI");
    var water = EVI.lt(0.1)
                   .and(MNDWI.gt(NDVI)//MNDWI>NDVI
                   .or(MNDWI.gt(EVI)));
    return image.addBands(water.rename("water"));
  }
  
  
  //Landsat-5/7处理方法
  var Landsat57 = {
    /*Landsat SR数据需要缩放，比例是 0.0001*/
    scaleImage: function(image) {
      var time_start = image.get("system:time_start");
      image = image.select(["B1","B2","B3","B4","B5","B7"]);
      image = image.divide(10000);
      image = image.set("system:time_start", time_start);
      return image;
    },
  
    /* SR数据去云*/
    srCloudMask: function(image) {
      var qa = image.select('pixel_qa');
      var cloudShadowBitMask = (1 << 3);
      var snowBitMask = (1 << 4);
      var cloudsBitMask = (1 << 5);
      var mask1 = qa.bitwiseAnd(cloudsBitMask).eq(0)
                    .and(qa.bitwiseAnd(snowBitMask).eq(0))
                    .and(qa.bitwiseAnd(cloudShadowBitMask).eq(0));
      var mask2 = image.mask().reduce(ee.Reducer.min());
      return image.updateMask(mask1.and(mask2));
    },
  
    //NDVI
    NDVI: function(image) {
        return image.addBands(image.normalizedDifference(["B4", "B3"])
                                   .rename("NDVI"));
    },
  
    //MNDWI
    MNDWI: function(image) {
        return image.addBands(image.normalizedDifference(["B2", "B5"])
                                   .rename("MNDWI"));
    },
  
    // EVI
    EVI: function(image) {
      var evi = image.expression("EVI = 2.5 * (NIR - R) / (NIR + 6*R -7.5*B + 1)", {
        NIR: image.select("B4"),
        R: image.select("B3"),
        B: image.select("B1")
      });
      return image.addBands(evi);
    },
  
    /*获取Landsat5的SR数据*/
    getL5SRCollection : function(startDate, endDate, region) {
      var dataset = l5_sr.filterDate(startDate, endDate)
                        .filterBounds(region)
                        .map(Landsat57.srCloudMask)
                        .map(Landsat57.scaleImage)
                        .map(Landsat57.NDVI)
                        .map(Landsat57.MNDWI)
                        .map(Landsat57.EVI)
                        .map(calcWater)
                        .select("water");
      return dataset;
    },
  
    /** 获取Landsat7的SR数据*/
    getL7SRCollection : function(startDate, endDate, region) {
      var dataset = l7_sr.filterDate(startDate, endDate)
                        .filterBounds(region)
                        .map(Landsat57.srCloudMask)
                        .map(Landsat57.scaleImage)
                        .map(Landsat57.NDVI)
                        .map(Landsat57.MNDWI)
                        .map(Landsat57.EVI)
                        .map(calcWater)
                        .select("water");
      return dataset;
    }
  };
  
  //Landsat 8处理方法
  var Landsat8 = {
    /* Landsat8 SR数据缩放*/
    scaleImage: function(image) {
      var time_start = image.get("system:time_start");
      image = image.select(["B2","B3","B4","B5","B6","B7"]);
      image = image.divide(10000);
      image = image.set("system:time_start", time_start);
      return image;
    },
  
    /*SR数据去云*/
    srCloudMask: function(image) {
      var cloudShadowBitMask = (1 << 3);
      var snowBitMask = (1 << 4);
      var cloudsBitMask = (1 << 5);
      var qa = image.select('pixel_qa');
      var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                   .and(qa.bitwiseAnd(snowBitMask).eq(0))
                   .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
      return image.updateMask(mask);
    },
  
    //NDVI
    NDVI: function(image) {
        return image.addBands(image.normalizedDifference(["B5", "B4"])
                                   .rename("NDVI"));
    },
  
    //MNDWI
    MNDWI: function(image) {
        return image.addBands(image.normalizedDifference(["B3", "B6"])
                                   .rename("MNDWI"));
    },
  
    // EVI
    EVI: function(image) {
      var evi = image.expression("EVI = 2.5 * (NIR - R) / (NIR + 6*R -7.5*B + 1)", {
        NIR: image.select("B5"),
        R: image.select("B4"),
        B: image.select("B2")
      });
      return image.addBands(evi);
    },
  
    /*获取Landsat8的SR数据*/
    getL8SRCollection : function(startDate, endDate, roi) {
      var dataset = l8_sr.filterDate(startDate, endDate)
                        .filterBounds(roi)
                        .map(Landsat8.srCloudMask)
                        .map(Landsat8.scaleImage)
                        .map(Landsat8.NDVI)
                        .map(Landsat8.MNDWI)
                        .map(Landsat8.EVI)
                        .map(calcWater)
                        .select("water");
      return dataset;
    }
  };
  
  //export
  function exportImageToDrive(image, key, region) {
    Export.image.toDrive({
      image: image, 
      description: "Drive-"+key,
      fileNamePrefix: key, 
      region: region,
      scale: 30,
      maxPixels: 1e13
    });
  }
  
  //去掉阴影
  function removeShadow(image, region) {
    var hand = ee.ImageCollection('users/gena/global-hand/hand-100');
    var hand30 = hand.mosaic().focal_mean(0.1).rename('elevation');
    var hillShadowMask = hand30.select('elevation').lte(50);
    var waterMask = image.updateMask(hillShadowMask.and(image.gte(0.8)))
                         .gte(0.25)
                         .clip(region);
    waterMask = waterMask.connectedPixelCount(50, true);
    waterMask = waterMask.updateMask(waterMask.gte(50));
    return image.updateMask(waterMask);
  }
  
  //生成每一年的水体 and输出像元值 
  // 输入参数为年份和区域
  function processYearWaterImage(year, month, region) {

  var startDate = ee.Date.fromYMD(year, month, 1);  
  var nextMonth = ee.Number(month).add(1);
  var endDate = ee.Date.fromYMD(year, nextMonth, 1); 
  
    var l5Water = Landsat57.getL5SRCollection(startDate, endDate, region);
    var l7Water = Landsat57.getL7SRCollection(startDate, endDate, region);
    var l8Water = Landsat8.getL8SRCollection(startDate, endDate, region);
    var waterImgs = l5Water.merge(l7Water).merge(l8Water)
    
    /*计算水体的频率*/
    var waterImg = waterImgs.sum()
                            .divide(waterImgs.count())
                            .clip(region);
    /*大于0.25像素为水体*/
    var hand = ee.ImageCollection('users/gena/global-hand/hand-100');
    var hand30 = hand.mosaic().focal_mean(0.1).rename('elevation');
    var hillShadowMask = hand30.select('elevation').lte(50);
    waterImg = waterImg.updateMask(hillShadowMask).updateMask(waterImg.gte(0.25));//mask外 is NoData
    waterImg = removeShadow(waterImg,roi_1);
    waterImg = waterImg.unmask(0).clip(region);//mask外is 0
    waterImg = waterImg.where(waterImg.select([0]).gt(0).and(waterImg.select([0]).lte(0.45)),0.45)
                      .where(waterImg.select([0]).gt(0.45).and(waterImg.select([0]).lte(0.5)),0.5)
                      .where(waterImg.select([0]).gt(0.5).and(waterImg.select([0]).lte(0.55)),0.55)
                      .where(waterImg.select([0]).gt(0.55),1);
  
    var key = "landsatWater-"+year+"_"+month;
    Map.addLayer(waterImg, {min:0,max:1,palette:['000000','blue']}, "water"+key, true);
    exportImageToDrive(waterImg, key, region);
    print('waterImgs_'+year+"_"+month,waterImgs);
      var stats2 = waterImg.reduceRegion({
        reducer: ee.Reducer.sum(),
        geometry:roi_1,
        scale: 30,
        maxPixels: 1E13
      });
    print(stats2,year)
  
  }
  
  //循环导出所有的水体
  function main() {
    //Map.addLayer(roi.style({color: "ffff00", fillColor: "00000000"}), {}, "roi");
    //开始年份和结束年份
    var startYear = 2019;
    var endYear = 2019;
    for (var year=startYear; year<=endYear; year++) {
      for (var month = 4 ; month<7;month++){
        processYearWaterImage(year,month ,roi_1)
      }
    }
  }
  
  main();