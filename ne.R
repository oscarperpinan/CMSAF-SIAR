library(maptools)

neAdm <- readShapePoly('~/Datos/NaturalEarth/ne_50m_admin_0_countries',
                       proj4string=CRS('+proj=longlat +ellps=WGS84'))

##spplot(neAdm["level"])

## France, Andorra, Portugal and Morocco
neighbours <- neAdm[neAdm$name %in% c('France', 'Andorra', 'Portugal', 'Morocco'), ]

writePolyShape(neighbours, 'data/neighbours.shp')
