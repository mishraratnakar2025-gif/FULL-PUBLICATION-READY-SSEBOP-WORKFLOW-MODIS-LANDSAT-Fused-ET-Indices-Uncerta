# FULL-PUBLICATION-READY-SSEBOP-WORKFLOW-MODIS-LANDSAT-Fused-ET-Indices-Uncerta
FULL PUBLICATION-READY SSEBOP WORKFLOW   MODIS + LANDSAT + Fused ET   Indices, Uncertainty, Confidence, Quadrants, Charts



/********************************************************
  FULL PUBLICATION-READY SSEBOP WORKFLOW
  MODIS + LANDSAT + Fused ET
  Indices, Uncertainty, Confidence, Quadrants, Charts
********************************************************/

/* ================================
   1️⃣ ROI
================================ */
var roi = ee.FeatureCollection('projects/eng-charge-478806-c6/assets/GOPAD_BASIN');

/* ================================
   2️⃣ SAFE HELPER FUNCTIONS
================================ */
function safeMean(ic, band){
  var img = ic.select(band).reduce(ee.Reducer.mean());
  return ee.Image(ee.Algorithms.If(img.bandNames().size().gt(0), img.rename(band), ee.Image.constant(0).rename(band)));
}
function safeStd(ic, band){
  var img = ic.select(band).reduce(ee.Reducer.stdDev());
  return ee.Image(ee.Algorithms.If(img.bandNames().size().gt(0), img.rename(band), ee.Image.constant(0).rename(band)));
}

/* ================================
   3️⃣ MODIS ET 2001–2023
================================ */
var modisET = ee.ImageCollection('MODIS/061/MOD16A2')
  .filterBounds(roi)
  .filterDate('2001-01-01','2023-12-31')
  .select('ET')
  .map(function(img){
    return img.multiply(0.1)
      .rename('ET_Monsoon')
      .copyProperties(img, ['system:time_start']);
  });

function addTime(img){
  var d = img.date();
  return img.set({
    'year': d.get('year'),
    'month': d.get('month')
  });
}
var modisET_time = modisET.map(addTime);

var years = ee.List.sequence(2001,2023);

// Function to compute seasonal ET
function seasonalET(monthList, name){
  return ee.ImageCollection(
    years.map(function(y){
      var ic = modisET_time.filter(ee.Filter.eq('year', y))
                            .filter(ee.Filter.inList('month', monthList));
      var img = ee.Image(ee.Algorithms.If(ic.size().gt(0),
                                         ic.sum().rename(name).toFloat(),
                                         ee.Image.constant(0).rename(name).toFloat()));
      return img.set('year', y)
                .set('system:time_start', ee.Date.fromYMD(y,6,1).millis());
    })
  );
}

// Seasonal ET
var ET_preMon  = seasonalET([3,4,5],   'ET_PreMonsoon');
var ET_monsoon = seasonalET([6,7,8,9], 'ET_Monsoon');
var ET_postMon = seasonalET([10,11],   'ET_PostMonsoon');
var ET_winter  = seasonalET([12,1,2],  'ET_Winter');

// Display
Map.addLayer(safeMean(ET_monsoon,'ET_Monsoon').clip(roi),
  {min:300,max:900,palette:['white','cyan','blue']},
  'Monsoon ET');

/* ================================
   4️⃣ CHIRPS SPI-12
================================ */
var rain = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
  .filterBounds(roi)
  .filterDate('2001-01-01','2023-12-31')
  .select('precipitation')
  .map(function(img){ return img.rename('P'); });

var rain_time = rain.map(addTime);

// Compute Monsoon rainfall per year (JJAS)
var monsoonRain = ee.ImageCollection(
  years.map(function(y){
    var ic = rain_time.filter(ee.Filter.eq('year', y))
                      .filter(ee.Filter.inList('month', [6,7,8,9]));
    var img = ee.Image(ee.Algorithms.If(ic.size().gt(0),
                                       ic.sum().rename('P'),
                                       ee.Image.constant(0).rename('P')));
    return img.set('year',y);
  })
);

// Monsoon SPI
var monRainMean = safeMean(monsoonRain,'P');
var monRainStd  = safeStd(monsoonRain,'P').max(0.001);

var SPI_monsoon = monsoonRain.map(function(img){
  return img.subtract(monRainMean)
            .divide(monRainStd)
            .rename('SPI_monsoon')
            .set('year', img.get('year'));
});

/* ================================
   5️⃣ MONSOON ET TREND (SEN SLOPE)
================================ */
// Add time band for Sen slope
var ET_monsoon_t = ET_monsoon.map(function(img){
  var year = ee.Number(img.get('year'));
  return img.addBands(ee.Image.constant(year).rename('t').toFloat());
});

// Compute Sen slope
var ET_monsoon_sen = ET_monsoon_t
  .select(['t','ET_Monsoon'])
  .reduce(ee.Reducer.sensSlope())
  .select('slope')
  .rename('ET_Monsoon_SenSlope')
  .clip(roi);

Map.addLayer(
  ET_monsoon_sen,
  {min:-10, max:10, palette:['red','white','blue']},
  'Monsoon ET Trend (Sen slope)'
);

// Optional significance mask
var sigMask = ET_monsoon_sen.abs().gt(5);  // mm/year threshold
Map.addLayer(
  ET_monsoon_sen.updateMask(sigMask),
  {min:-10, max:10, palette:['red','white','blue']},
  'Monsoon ET Significant Trend'
);

/* ================================
   6️⃣ MONSOON DROUGHT YEARS (SPI-12)
================================ */
var SPI12_annual = ee.ImageCollection(
  years.map(function(y){
    var ic = SPI_monsoon.filter(ee.Filter.eq('year', y));
    return safeMean(ic,'SPI_monsoon').set('year',y);
  })
);

var droughtYears = SPI12_annual.filter(ee.Filter.lt('SPI_monsoon', -1))
                               .aggregate_array('year');
var normalYears  = SPI12_annual.filter(ee.Filter.gt('SPI_monsoon', -0.5))
                               .aggregate_array('year');

/* ================================
   7️⃣ MONSOON RESISTANCE–RESILIENCE
================================ */
var ET_monsoon_drought = ET_monsoon.filter(ee.Filter.inList('year', droughtYears));
var ET_monsoon_normal  = ET_monsoon.filter(ee.Filter.inList('year', normalYears));

var mean_drought = safeMean(ET_monsoon_drought,'ET_Monsoon');
var mean_normal  = safeMean(ET_monsoon_normal,'ET_Monsoon').max(0.001);
var std_drought  = safeStd(ET_monsoon_drought,'ET_Monsoon').max(0.001);

var Resistance = mean_drought.divide(mean_normal)
  .rename('Resistance')
  .clip(roi);
var Resilience = mean_normal.divide(std_drought)
  .rename('Resilience')
  .clip(roi);

Map.addLayer(Resistance,
  {min:0,max:1,palette:['red','yellow','green']},
  'Monsoon ET Resistance');

Map.addLayer(Resilience,
  {min:0,max:3,palette:['red','yellow','green']},
  'Monsoon ET Resilience');

/* ================================
   8️⃣ MODIS–LANDSAT FUSION 2019–2023
================================ */
var landsat = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(roi)
  .filterDate('2019-01-01','2023-12-31')
  .map(function(img){
    var qa = img.select('QA_PIXEL');
    var mask = qa.bitwiseAnd(1<<3).eq(0).and(qa.bitwiseAnd(1<<4).eq(0));
    var LST = img.select('ST_B10').multiply(0.00341802).add(149.0).rename('LST');
    return LST.updateMask(mask).copyProperties(img,['system:time_start']);
  });

var LSTmin = landsat.reduce(ee.Reducer.percentile([5])).rename('LSTmin');
var LSTmax = landsat.reduce(ee.Reducer.percentile([95])).rename('LSTmax');

var EF = landsat.map(function(img){
  return LSTmax.subtract(img).divide(LSTmax.subtract(LSTmin)).clamp(0,1).rename('EF')
    .copyProperties(img,['system:time_start']);
});

var Rn = 15;   // MJ m-2 day-1
var lambda = 2.45; // MJ kg-1
var ET_Landsat = EF.map(function(img){ return img.multiply(Rn).divide(lambda).rename('ET').copyProperties(img,['system:time_start']); });

var annualET_LS = ee.ImageCollection(
  ee.List.sequence(2019,2023).map(function(y){
    var ic = ET_Landsat.filter(ee.Filter.calendarRange(y,y,'year'));
    return safeMean(ic,'ET').multiply(365).rename('ET').set('year',y);
  })
);

// 🔹 RENAME MODIS bands to ET for fusion
var MODIS_recent = ET_monsoon.filter(ee.Filter.inList('year',[2019,2020,2021,2022,2023]))
                             .map(function(img){ return img.rename('ET'); });

var fusedET = ee.ImageCollection(
  ee.List.sequence(2019,2023).map(function(y){
    var mod = MODIS_recent.filter(ee.Filter.eq('year',y)).first();
    var ls  = annualET_LS.filter(ee.Filter.eq('year',y)).first();
    return ee.Image(mod).multiply(0.6)
           .add(ee.Image(ls).multiply(0.4))
           .rename('ET')
           .set('year',y);
  })
);

var fusedET_mean = safeMean(fusedET,'ET').clip(roi);
Map.addLayer(fusedET_mean,{min:500,max:1200,palette:['white','cyan','blue']},'Fused MODIS–Landsat ET');

/* ================================
   9️⃣ EXPORTS
================================ */
var exportList = [
  {img:Resistance, name:'ET_Resistance_500m'},
  {img:Resilience, name:'ET_Resilience_500m'},
  {img:fusedET_mean, name:'Fused_ET_MODIS_Landsat_250m'},
  {img:ET_monsoon_sen, name:'Monsoon_ET_SenSlope_500m'}
];

exportList.forEach(function(o){
  Export.image.toDrive({
    image: o.img,
    description: o.name,
    region: roi,
    scale: 500,
    crs:'EPSG:4326',
    maxPixels:1e13
  });
});

/* ================================
   🔟 ANNUAL ET TIME SERIES CHART
================================ */
function setTimeStart(ic){ 
  return ic.map(function(img){ 
    return img.set('system:time_start', ee.Date.fromYMD(ee.Number(img.get('year')),1,1).millis()); 
  }); 
}

var MODIS_ts = setTimeStart(MODIS_recent);
var LS_ts    = setTimeStart(annualET_LS);
var Fused_ts = setTimeStart(fusedET);

var combinedChart = ui.Chart.image.seriesByRegion({
  imageCollection: MODIS_ts.merge(LS_ts).merge(Fused_ts),
  band: 'ET',
  regions: roi,
  reducer: ee.Reducer.mean(),
  scale: 500,
  seriesProperty: 'year'
}).setOptions({
  title:'Annual ET Comparison (MODIS, Landsat, Fused)',
  hAxis:{title:'Year'},
  vAxis:{title:'ET (mm/year)'},
  lineWidth:2,
  pointSize:4,
  colors:['blue','green','red']
});
print(combinedChart);



////////////////

/* ================================
   11️⃣ MONSOON RESISTANCE–RESILIENCE QUADRANTS
================================ */

// Compute mean values across ROI
var Rm = Resistance.reduceRegion({
  reducer: ee.Reducer.mean(), 
  geometry: roi, 
  scale: 500, 
  maxPixels: 1e13
}).getNumber('Resistance');

var Sm = Resilience.reduceRegion({
  reducer: ee.Reducer.mean(), 
  geometry: roi, 
  scale: 500, 
  maxPixels: 1e13
}).getNumber('Resilience');

// Quadrants:
// 0: Low Resistance, Low Resilience
// 1: Low Resistance, High Resilience
// 2: High Resistance, Low Resilience
// 3: High Resistance, High Resilience
var quadrant = Resistance.gt(Rm).multiply(2)
                .add(Resilience.gt(Sm))
                .rename('Quadrant')
                .clip(roi);

Map.addLayer(
  quadrant,
  {min:0, max:3, palette:['red','yellow','blue','darkgreen']},
  'Resistance–Resilience Quadrant'
);

/* ================================
   12️⃣ LEGEND FUNCTION
================================ */
function addLegend(title, palette, min, max) {
  var legend = ui.Panel({style:{position:'bottom-left', padding:'8px 15px', backgroundColor:'white'}});
  var legendTitle = ui.Label(title, {fontWeight:'bold', fontSize:'14px', margin:'0 0 4px 0'});
  legend.add(legendTitle);
  var lon = ee.Image.pixelLonLat().select('latitude');
  var gradient = lon.multiply((max-min)/100.0).add(min);
  var colorBar = ui.Thumbnail({
    image: gradient.visualize({min:min, max:max, palette:palette}), 
    style:{stretch:'horizontal', height:'12px', margin:'0 0 4px 0'}
  });
  legend.add(colorBar);
  var panel = ui.Panel({layout: ui.Panel.Layout.flow('horizontal')});
  panel.add(ui.Label(min));
  panel.add(ui.Label((min+max)/2, {margin: '0 45px'}));
  panel.add(ui.Label(max));
  legend.add(panel);
  Map.add(legend);
}

// Add quadrant legend
addLegend('Resistance–Resilience Quadrant', ['red','yellow','blue','darkgreen'], 0, 3);

/* ================================
   13️⃣ EXPORT QUADRANT MAP
================================ */
Export.image.toDrive({
  image: quadrant,
  description: 'ET_Quadrants_500m',
  region: roi,
  scale: 500,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});
