using DataFrames
using Gadfly
svg = true

tags = [
"CPAR100000lin49x50T49.0TC", 
"CPAR100000lin49x20T49.0TC",
"CPAR100000lin49x10T49.0TC",
]

function rmean(v)
    cumsum(v)./[1:length(v)]
end

theme = Theme(panel_fill=color("light grey"), grid_color=color("white"))
geo = Geom.line

for tag in tags
    theme = Theme(panel_fill=color("light grey"), grid_color=color("white"), default_color=color("black"),
    default_point_size=0.01cm)
    geo = Geom.line

    thetas = readcsv("$tag/thetas$tag.csv")'
    N = min(size(thetas,2), 100000)
    thetas = thetas[:, 1:10:N]
    N = size(thetas,2)


    dfth = DataFrame(n=1:N, c1=thetas[1,:][:], c2=thetas[2,:][:], c3=thetas[3,:][:], c4=thetas[4,:][:], 
                            c5=thetas[5,:][:], c6=thetas[6,:][:], c7=thetas[7,:][:], c8=thetas[8,:][:], 
    tr="iterates")
    dftrue = DataFrame(n=1:N, c1=0.1, c2=0.7, c3=0.35, c4=0.2, c5=0.1, c6=0.9, c7=0.3, c8=0.1, tr="Truth")
    dfrmean = DataFrame(n=1:N, c1=rmean(thetas[1,:][:]), c2=rmean(thetas[2,:][:]), c3=rmean(thetas[3,:][:]), c4=rmean(thetas[4,:][:]), 
                               c5=rmean(thetas[5,:][:]), c6=rmean(thetas[6,:][:]), c7=rmean(thetas[7,:][:]), c8=rmean(thetas[8,:][:]), 
    tr="Running mean")
    
    
    dftrue2 = stack(dftrue, [2,3,4,5,6,7,8,9])
    dfth2 = stack(dfth, [2,3,4,5,6,7,8,9])
    dfrm2 = stack(dfrmean, [2,3,4,5,6,7,8,9])
    df = vcat(dftrue2, dfth2, dfrm2)
 
    p = plot(df, x="n", y="value", color="tr", ygroup="variable", Geom.subplot_grid(Stat.yticks,  geo, free_y_axis=true),
    Guide.xlabel("Iteration"), Guide.ylabel("Parameters"), Guide.title("Trace"),  Guide.colorkey(""),
    Scale.discrete_color_manual("green","red","grey12"), Scale.x_continuous(format=:plain), Scale.y_continuous(format=:plain),
    theme
    )
    draw(PNG("$tag/thetas$tag.png",18cm, 22cm), p) 
    draw(PDF("$tag/thetas$tag.pdf",18cm, 22cm), p)
    

end


