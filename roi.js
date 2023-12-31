var input = ee.Image("JRC/GSW1_4/GlobalSurfaceWater"),
    roi = 
    /* color: #ffc82d */
    /* shown: false */
    ee.Geometry.Polygon(
        [[[122.5521613490698, 39.971425268402754],
          [122.54357828022215, 39.965768158244494],
          [122.5372264789685, 39.958991703924546],
          [122.53602517963621, 39.95050478198539],
          [122.51491083027098, 39.96813630783827],
          [122.50323785663817, 39.961821059956264],
          [122.50598443866942, 39.9481360217183],
          [122.49568475605223, 39.9481360217183],
          [122.4888183009741, 39.93655423522644],
          [122.50186456562254, 39.91865125475899],
          [122.47851861835692, 39.93076260683187],
          [122.46959222675535, 39.92075773075954],
          [122.48675836445067, 39.90548431169918],
          [122.47439874531004, 39.90390410844596],
          [122.4888183009741, 39.88862693140338],
          [122.47165216327879, 39.88151398047327],
          [122.46547235370848, 39.879142832845375],
          [122.45843427785928, 39.877694392502484],
          [122.45173960579609, 39.87624549844747],
          [122.44178308368895, 39.87650812820864],
          [122.44830621601317, 39.86649532746978],
          [122.47439874531004, 39.86438718318194],
          [122.48195184589598, 39.87176540478728],
          [122.49293817402098, 39.84224775745268],
          [122.48195184589598, 39.84013886828907],
          [122.47302545429442, 39.858062361663706],
          [122.4448729884741, 39.86333308674682],
          [122.44075311542723, 39.854372613227326],
          [122.45997918964598, 39.847519697003115],
          [122.4555159938452, 39.83829353714216],
          [122.43457330585692, 39.844883777828166],
          [122.42702020527098, 39.839611635878356],
          [122.40161432148192, 39.8359208956635],
          [122.38376153827879, 39.825902172432265],
          [122.41054071308348, 39.81640940259391],
          [122.44658960224363, 39.813772289360806],
          [122.43697656513426, 39.8026953089109],
          [122.44933618427488, 39.79636480449595],
          [122.47474206806395, 39.8063878340335],
          [122.47611535907957, 39.815618279248014],
          [122.51353753925535, 39.82194701103332],
          [122.5327636134741, 39.83381181220331],
          [122.52383722187254, 39.85147338601085],
          [122.52280725361078, 39.856480820271834],
          [122.5217772853491, 39.863860136988635],
          [122.5437499415991, 39.860697774810596],
          [122.53482354999754, 39.87334634903327],
          [122.5547362697241, 39.86913041673659],
          [122.5547362697241, 39.88493882674066],
          [122.53345025898192, 39.8891537873016],
          [122.53070367695067, 39.898636501051115],
          [122.5712157619116, 39.89600254550048],
          [122.59524835468504, 39.91022470278069],
          [122.5931884181616, 39.92444390783663],
          [122.57739557148192, 39.92023111783429],
          [122.56640924335692, 39.92655020563699],
          [122.55061639667723, 39.92391732326138],
          [122.53551019550535, 39.92497048836182],
          [122.57876886249754, 39.93497474881329],
          [122.57327569843504, 39.94760961940852],
          [122.55816949726317, 39.957610570620695],
          [122.5602294337866, 39.968662552162925]]]);
var visualization = {
  bands: ['max_extent'],
  min: 0,
  max: 1,
  palette: ['green', 'blue']
};
Map.setCenter(122.51, 39.88, 8);
Map.addLayer(input, visualization, 'max_extent');
Map.addLayer(roi, {}, 'max_extent');
var output = input.select('max_extent').clip(roi);
output = output.rename(['max_extent']);
print(output,"output")
Map.addLayer(output, visualization, 'output');

// 将 "max_extent" 波段设为二值图像，其中值为 1 的像素为 1，其他像素为 0
var binaryImage = output.select('max_extent').eq(1);
// 计算值为 1 的像素的总面积
var area = binaryImage.multiply(ee.Image.pixelArea()).reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: roi,
  scale: 30,
  maxPixels: 1e13
});
// 通过 getInfo 方法获取字典的 JavaScript 对象
var areaValue = area.get('max_extent');
// 打印结果
print('Total area with value 1:', ee.Number(areaValue).divide(1e6)); // 将平方米转换为平方千米

//导出roi
var feature = ee.Feature(roi);
var featureCollection = ee.FeatureCollection([feature]);

// 导出为 Shapefile
Export.table.toDrive({
  collection: featureCollection,
  description: 'drawn_polygon_export',
  folder: 'GEE_exports',
  fileNamePrefix: 'BLHreservoir',
  fileFormat: 'SHP'
});
Export.image.toDrive({
  image: output,
  description: 'maxWaterExtent',
  folder: 'GEE_exports',
  fileNamePrefix: 'maxbiliuheWaterExtent',
  scale: 30,
  region: roi,
  maxPixels: 1e13,
  crs: 'EPSG:4326',
});