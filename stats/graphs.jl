using CSV
using StatsPlots
using DataFrames
using Plots

df = DataFrame(CSV.File("data.csv", delim = '\t'));

# Divide all the values by 1000 to turn them into ms.
for i = 2:6
  df[!, i] = df[!, i] ./ 1000;
end
print(df);
print("\n\n");

names = vcat(["sequential"], repeat(["pthread", "openMP", "openCilk"], inner = 3));
groups = vcat(["1 thread"], repeat(["2 threads", "4 threads", "8 threads"], outer = 3));

plotlyjs()
b1 = groupedbar(
  names, df.belgium_osm, group = groups,
  title = "belgium_osm",
  xlabel = "Algorithm used",
  ylabel = "Execution times in ms"
);

b2 = groupedbar(
  names, df.dblp2010, group = groups,
  title = "dblp-2010",
  xlabel = "Algorithm used",
  ylabel = "Execution times in ms"
);

b3 = groupedbar(
  names, df.NACA0015, group = groups,
  title = "NACA0015",
  xlabel = "Algorithm used",
  ylabel = "Execution times in ms"
);

b4 = groupedbar(
  names, df.mycielskian13, group = groups,
  title = "mycieskian13",
  xlabel = "Algorithm used",
  ylabel = "Execution times in ms"
);

b5 = groupedbar(
  names, df.comYoutube, group = groups,
  title = "com-Youtube",
  xlabel = "Algorithm used",
  ylabel = "Execution times in ms"
);