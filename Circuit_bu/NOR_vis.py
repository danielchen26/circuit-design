# --------- ########## --------- #########-------------------   visualizations  --------- ########## --------- #########-------------------
# Import visualizations modules
import pandas as pd
import IPython
import seaborn as sns
import altair as alt
def vegify(spec):
    IPython.display.display({
        'application/vnd.vegalite.v2+json': spec.to_dict()
    }, raw=True)
alt.data_transformers.enable('default', max_rows=None)

# cd Circuit_bu  # change dir path

#  ------ First reimport data ------
DB_import = pd.read_csv('NOR_database/NOR_database_Corrected_withoutID.csv')
DB_import['fc'] = 1



DB_import.apply(lambda x: x)
# correct the optimial net variable names ------------------------------------
ext = DB_import['Expresso_tts'][222]
neti = DB_import['Optimal_net'][222]
expr.expr(ext).inputs


splt = [i+")" for i in neti.split("),")]
splt[-1] = splt[-1].replace("))", ")")
# splt
splt_modi = [j.replace("NOR", "Nor") for j in splt]
splt_array = np.array([i.split(" = ") for i in splt_modi])

net_dict = dict(zip(splt_array[:, 0], splt_array[:, 1]))
keys = list(net_dict.keys())
values = list(net_dict.values())

keys
values


r1 = net_dict['O']
r2 = r1.replace("I1", net_dict['I1'])


# get the final reduced exrepssion from network
output = values[0]
while any(x in output for x in keys):
    for ki in keys:
        if ki in output:
            print(ki)
            output = output.replace(ki, net_dict[ki])
            output
expression = expr.expr(output)
redu = expr2truthtable(expression)
# espresso_exprs(expr.expr(output))
tt = list(redu.iter_relation())
redu_br = ''.join([str(i[1]) for i in tt])

# correct the optimial net variable names ------------------------------------
DB_import[DB_import['Num_in'] == 3]
[expr.expr(i).inputs for i in DB_import['Expresso_tts']]
DB_import[DB_import['Num_in'] == 3]


# import plotly.graph_objs as go
#
# iplot([go.Histogram2dContour(x=DB_sum.Num_gates, y=DB_sum.fc, contours=dict(coloring='heatmap')),
#        go.Scatter(x=DB_sum.Num_gates, y=DB_sum.fc, mode='markers', marker=dict(color='white', size=3, opacity=0.3))], show_link=False)


p2 = alt.Chart(DB_import).mark_bar().encode(
    x='Num_gates',
    y='count(Num_gates)',
    color=alt.Color('Num_in:N', scale=alt.Scale(scheme='set1'))
).configure(background="White").interactive()
vegify(p2)


p3 = alt.Chart(DB_import.iloc[0:400, :]).mark_bar().encode(
    x='Num_in',
    y='count(Num_gates)',
    order=alt.Order('count(Num_gates)', sort='ascending'),
    color=alt.Color('Num_gates:N', scale=alt.Scale(scheme='set1'))
).configure(background="White").interactive()
vegify(p3)


alt.Chart(DB_import.iloc[0:400, :]).mark_bar(opacity=0.7).encode(
    x='Num_gates:O',
    y=alt.Y('count(Num_gates):Q', stack=None),
    color=alt.Color('Num_in:N', scale=alt.Scale(scheme='dark2')),
).configure(background="White")


# ---- normalized total number of functions with repect to each number of gate -----
# Relative abundance of different # of inputs for each min gates realization
relative_ab = alt.Chart(DB_import).mark_bar(opacity=0.7).encode(
    x=alt.X('Num_gates:O', title="# of min NOR gates"),
    y=alt.Y('count(Num_gates):Q', stack="normalize",
            title="Fraction of total # of boolean functions"),
    color=alt.Color('Num_in:N', scale=alt.Scale(
        scheme='dark2'), title='Number of inputs'),
).properties(title="Relative abundance of different # of inputs for each min gates realization"
             ).configure(background="White").interactive()

relative_ab.save('./Plots/relative_abundance.html')


alt.Chart(test).mark_bar(opacity=0.7).encode(
    x=alt.X('Num_in:O'),
    y=alt.Y('sum(fc):O', stack="normalize"),
    # order=alt.Order('sum(Num_gates)', sort='ascending'),
    color=alt.Color('Num_in:N', scale=alt.Scale(scheme='set1')),
    column='Num_gates:O',
).configure(background="White").interactive()


# --- Plot3 :  Total number of functions histogram for each min NOR net in different # input cases (for p-class: just change Y->cout())
hist1 = alt.Chart(DB_import).mark_bar().encode(
    x=alt.X('Num_gates:O', title="# of min NOR gates"),
    y=alt.Y('count(Num_gates):Q', title='Total # of boolean functions'),
    color=alt.Color('Num_gates:O', scale=alt.Scale(scheme='category20b')),
    column='Num_in:N').properties(title="Total number of functions histogram for each min NOR net in different # input cases"
                                  ).configure(background="White").interactive()
hist1.save('./Plots/hist_gates_in.html')


# --- Plot4 :  Total number of functions histogram(discrete dot: easy to read) for each min NOR net in different # input cases
hist2 = alt.Chart(DB_import).mark_bar().encode(
    x=alt.X('Num_gates:O', title="# of min NOR gates"),
    y=alt.Y('sum(fc):O', sort='ascending',
            title='Total # of boolean functions'),
    color=alt.Color('Num_gates:O', scale=alt.Scale(scheme='category20b')),
    column='Num_in:N'
).properties(title="Total number of functions square plot for each min NOR net in different # input cases"
             ).configure(background="White").interactive()
hist2.save('./Plots/hist_gates_in_square.html')


slider = alt.binding_range(min=1, max=14, step=1)
select_number_of_gates = alt.selection_single(
    name="NOR", fields=['Num_gates'], bind=slider)


slide_bar = alt.Chart(DB_import).mark_bar().encode(
    x=alt.X('Num_gates:N', axis=alt.Axis(title=None)),
    y=alt.Y('sum(fc):O', title='Total # of boolean functions'),
    color=alt.Color('Num_in:N', scale=alt.Scale(scheme='set1')),
    column='Num_in:O'
).properties(title="Total number of functions square plot for each min NOR net in different # input cases",
             width=150
             ).add_selection(
    select_number_of_gates
).transform_filter(
    select_number_of_gates
).configure(background="White")

slide_bar.save('./Plots/hist_gates_in_square_bar_control.html')


# -----------------------   interactive selection
from vega_datasets import data
# source = data.seattle_weather()
# scale = alt.Scale(domain=[2,3,4,5,6],
#                   range=['#e7ba52', '#a7a7a7', '#aec7e8', '#1f77b4', '#9467bd'])
# color = alt.Color('Num_gates:N', scale=scale)
color = alt.Color('Num_gates:O', scale=alt.Scale(
    scheme='tableau20'), legend=alt.Legend(title="Number gates"))
# We create two selections:
# - a brush that is active on the top panel
# - a multi-click that is active on the bottom panel
brush = alt.selection_interval(encodings=['x'])
click = alt.selection_multi(encodings=['color'])


norm_bar = alt.Chart().mark_bar(opacity=0.7).encode(
    x=alt.X('Num_gates:O', title='minimun NOR gates'),
    y=alt.Y('count(Num_gates):Q', stack="normalize",
            title='Total # of functions colored by inputs'),
    color=alt.Color('Num_in:O', legend=alt.Legend(title="Number inputs")),
).properties(
    width=600,
    height=300
).add_selection(
    brush
).transform_filter(
    click
).interactive()


# Bottom panel is a bar chart of weather type
bars = alt.Chart().mark_bar().encode(
    x=alt.X('Num_gates:N', title='minimun NOR gates'),
    y=alt.Y('count(Num_gates):Q', title='Total # of functions of all inputs'),
    color=alt.condition(click, color, alt.value('lightgray')),
).transform_filter(
    brush
).properties(
    width=600,
).add_selection(
    click
).interactive()

chart = alt.vconcat(
    norm_bar,
    bars,
    data=DB_import,
    title="NOR DATABASE"
).configure(background="White")

chart.save('./Plots/chart.html')


#  ---------  Interactive Crossfilter -------------
#

import altair as alt
from vega_datasets import data

source = alt.UrlData(
    data.flights_2k.url,
    format={'parse': {'date': 'date'}}
)

brush = alt.selection(type='interval', encodings=['x'])

# Define the base chart, with the common parts of the
# background and highlights
base = alt.Chart().mark_bar().encode(
    x=alt.X(alt.repeat('column'), type='quantitative'),
    y=alt.Y('count(Num_gates):Q', stack="normalize")
).properties(
    width=180,
    height=130
)

# blue background with selection
background = base.properties(selection=brush)

# yellow highlights on the transformed data
highlight = base.encode(
    color=alt.value('goldenrod')
).transform_filter(brush)

# layer the two charts & repeat
alt.layer(
    background,
    highlight,
    data=DB_import
).properties(
    width=1000,
    height=800
).repeat(column=["Num_in", "Num_gates"]).configure(background="White")


#  -------- Seaborn box and swarm plots -----------
import matplotlib.pyplot as plt
%config InlineBackend.figure_format = 'retina'

# style.use('ggplot')
plt.style.use('tableau-colorblind10')
sns.set(font_scale=1.6)
sns.set_style("white")
f1 = plt.figure(figsize=(12, 10))
sns.violinplot(x="Num_in", y="Num_gates", bw=0.1,
               scale="width", inner="quartile", data=DB_import)
# ax.xaxis.grid(True)
plt.xlabel("Number of inputs")
plt.ylabel("Number minimum NOR gates")
sns.despine(trim=True)
plt.tight_layout()
plt.axhline(y=7, color='r', linestyle=':')
plt.title(
    "Distribution of min NOR gates(Categorical[1~14]) of functions with different inputs[1~4] ")
sns.set_context("poster")
plt.show()

f1.savefig("./Plots/violin_each_inputs.png", dpi=600)
