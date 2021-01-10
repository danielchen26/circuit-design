using Plots, LaTeXStrings;gr()
#font = Plots.font("Helvetica", 12)
x = collect(0:0.1:10)
plt = plot(x, Hill_shift.(2.,0.2,2,5.,x),lw=3, color = :darkgreen, labels = L"Hill(x) = dn + \frac{(up - dn) * K^{n}}{(K^{n} + x^{n})}",
    ylims = (0,2.2), tickfont = Plots.font("Helvetica", 12), guidefont = Plots.font("Helvetica", 12), legendfontsize = 13, legend = :right, fg_legend = :transparent)

xlabel!("Input");ylabel!("Output")
scatter!(x, fill(0.2, size(x)), alpha = 0.55, color =:blue, markersize = 2,label ="dn")
scatter!(x, fill(2, size(x)), alpha = 0.55,color =:darkorange, markersize = 2,label ="up")
#annotate!((7.5, x[5], Plots.text("dn", 10, :red, :center)))
#savefig(plt, "SavedPlots/hill.png", dpi=300)