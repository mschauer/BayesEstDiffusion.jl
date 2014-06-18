using Gadfly
using DataFrames

function rmean(v)
cumsum(v)./[1:length(v)]
end
    
theme = Theme(panel_fill=color("light grey"), grid_color=color("white"))
geo = Geom.line
    
xobsdf = DataFrame(hcat(0:49,readcsv("autoreg50fo.csv"; has_header=false)))
names!(xobsdf, [:t,:RNA, :P, :P2, :DNA])
p = plot(melt(xobsdf, "t"), x="t", y="value", color="variable",
    Guide.ylabel("#"), Guide.title(""),
    theme )
draw(PNG("autoreg50fo.png",18cm, 9cm), p) 
draw(PDF("autoreg50fo.pdf",18cm, 9cm), p) 
             
             
             