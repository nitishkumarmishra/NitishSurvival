library(ggtern)
data('Feldspar')
ggplot(data=Feldspar,aes(Ab,An)) + 
  stat_density_2d(
    geom='raster',
    aes(fill=..density..),
    bins=5,
    contour = FALSE) +
  theme_bw()+
  #scale_fill_distiller(palette= "Spectral", direction=1) +
  scale_fill_distiller(palette= "RdBu", direction=1) +
  geom_point()
