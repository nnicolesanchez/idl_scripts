pro bargasfraction,plotter=plotter

loadct,4

if keyword_set(plotter) then paperplot,filename='bargasfraction_h277.ps'

;#### h277 [Galaxy, BH] Fractions, Respectively ####
earlyfrac	=[.07,.00]
shockedfrac	=[.13,.21]
coldfrac	=[.52,.48]
clumpyfrac	=[.28,.31]

bar_plot,earlyfrac+shockedfrac+coldfrac+clumpyfrac,colors=[220,220],title='h277 Gas Fractions',xtitle='Galaxy                  SMBH',ytitle='Gas Fraction'
bar_plot,shockedfrac+coldfrac+clumpyfrac,colors=[160,160],/overplot
bar_plot,coldfrac+clumpyfrac,colors=[50,50],/overplot
bar_plot,clumpyfrac,colors=[100,100],/overplot				;blue
xyouts,[.5,.5],'hello',/normal

if keyword_set(plotter) then paperplot,/close
stop
if keyword_set(plotter) then paperplot,filename='bargasfraction_h258.ps'

;#### h258 [Galaxy, BH] Fractions, Respectively ####
earlyfrac	=[.05,.02]
shockedfrac     =[.09,.09]
coldfrac        =[.48,.25]
clumpyfrac	=[.38,.64]

bar_plot,earlyfrac+shockedfrac+coldfrac+clumpyfrac,colors=[220,220],title='h258 Gas Fractions',xtitle='Galaxy                  SMBH',ytitle='Gas Fraction'
bar_plot,shockedfrac+coldfrac+clumpyfrac,colors=[160,160],/overplot
bar_plot,coldfrac+clumpyfrac,colors=[50,50],/overplot
bar_plot,clumpyfrac,colors=[100,100],/overplot

if keyword_set(plotter) then paperplot,/close


end

