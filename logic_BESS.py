import pandas as pd
import json
from pandas.io.json import json_normalize
import itertools
import numpy as np
from pyeda.inter import *
from graphviz import Source
import re, pyeda.boolalg.expr as expr
import logic_fun as lf

# ---------------------------  Importing BESS NAND database (inputs : 2~4 ) -----------------------------
# library of optimial gates
df2 = pd.read_csv("database/2NAND.csv")
df2 = df2.replace(r'\r\n', ',', regex = True)# modify df[['Optimal Network']] column
df3 = pd.read_csv("database/3NAND.csv")
df3 = df3.replace(r'\r\n', ',', regex = True)# modify df[['Optimal Network']] column
df4 = pd.read_csv("database/BESS_NAND.csv")
df4 = df4.replace(r'\r\n', ',', regex = True)# modify df[['Optimal Network']] column


# ------------------- Generating a dictionary of minimal dnf design of each "Num_fun" of boolean functions
Num_fun, Boolean_dict = lf.Boolean_expresso_gen()
# convertion dict which switches 'Or' and 'And'
convertion = dict([('Or', 'And'),('And', 'Or')])

# Creating NOR optimal net database(takes quite a bit of time!!!)
database = []
t = 1
for func, expression in Boolean_dict.items():
    print("The {}\n boolean function {}\n has minimal tts: {}".format(t, func, expression) )

    # converted expression
    func_converted = expr.expr(lf.multiple_replace(convertion, str(expression[0])))

    # Final output of the NOR gates optimal network of the input boolean fucntion (1."And" "Or" Switched )
    Opt_NOR_net = lf.permu_search_gen(func_converted,df4)

    data = np.append([func,expression],Opt_NOR_net)
    database.append(list(data))
    t+=1
    if t>10:break


NOR_database = np.array(database)

NOR_df = pd.DataFrame({'Binary':NOR_database[:,0],'Expresso_tts':np.array([i[0] for i in NOR_database[:,1]]),'Num_gates':NOR_database[:,2],'Optimal_net':NOR_database[:,3]})
# NOR_df.to_csv("NOR_df.csv", index=False)



#
# # ---------------------------------------------------------------------
#
# # Test specific boolean function
# st = '1111111100000000'
# st = '0000001101010111'
# express = Boolean_dict[st]# ~d
# func_converted = expr.expr(lf.multiple_replace(convertion, str(express[0])))
#
# import pyeda.boolalg.expr as expr
# all_permu = list(itertools.permutations('abcd',4))
# origin = ('a','b','c','d')
# pm_dict_list = [dict(zip(origin, each)) for each in all_permu]
#
# # assinging input boolean function
# rs_result = func_converted
#
# # The following block test each permutation case in the permutation dictionary list,
# # totally 24 possible case, and choose one the match integer that corresponds to the data library
#
# NOR = []
# for per_num in range(24):
#     #
#     rs_exchange = lf.multiple_replace(pm_dict_list[per_num], str(rs_result))
#     rs_exchange
#
#     # 3 possible correction terms, only one is needed to be corrected.
#     correction = {'Anc': 'And', 'Anb': 'And', 'Ana': 'And'}
#     for key in correction.keys():
#         if key in rs_exchange:
#             print(key)
#             correction_sub = {key: correction[key]}
#             rs_ex_permu = expr.expr(lf.multiple_replace(correction_sub,rs_exchange))
#         else:
#             rs_ex_permu = expr.expr(lf.multiple_replace(correction,rs_exchange))
#
#
#     # rs_ex_permu =  expr.expr(multiple_replace(correction,rs_exchange))
#     rs_ex_permu_tt = list (rs_ex_permu.iter_relation())
#     permu_br = ''.join([str(i[1]) for i in rs_ex_permu_tt])
#     permu_intr = int(permu_br,2)
#
#     if not df[df['Integer'] == permu_intr].empty:
#             opt_net = df[df['Integer'] == permu_intr]
#             opt_struct = opt_net.values[0][3]
#             opt_num_gates = opt_net.values[0][2]
#             print('The original network structure is :\n', opt_struct)
#
#             print('we have use the permuation dictionary of ', pm_dict_list[per_num])
#
#             # generate reverse permuation mapping
#             pm_reverse = dict((v,k) for k, v in pm_dict_list[per_num].items())
#             print('The reversed permutation mapping is:\n', pm_reverse)
#
#             # construct the new NOR gate library
#             NOR_net = lf.multiple_replace(pm_reverse, opt_struct).replace('NAND','NOR')
#             print("The optimal NOR network structure is:\n", NOR_net)
#
#             NOR = [opt_num_gates, NOR_net]
#             break
#
#     print(df[df['Integer'] == permu_intr])
# # ---------------------------------------------------------------------
#
#
#
# # ------- Some visualizations ------
# import matplotlib.pyplot as plt
# import altair as alt
#
# NOR_reload =  pd.read_csv("NOR_df2.csv")
# type(NOR_reload.Binary[3])
# NOR_reload['Binary'] = list(Boolean_dict.keys())# somehow need to force the string type
# columns = NOR_reload.columns
#
# NOR_stat = NOR_reload.groupby('Num_gates').count().reset_index().replace('NAN',np.nan).dropna()
# NOR_stat['Num_gates'] = NOR_stat['Num_gates'].astype(int)
# NOR_stat.sort_values(by=['Num_gates'])#.apply(lambda x: x.cumsum())
# # NOR_reload.to_csv("NOR_df2.csv", index =False)
# NOR_stat
# NOR_reload[NOR_reload['Binary'] =='1111111100000000']






# -----------------------------------------------    Generating all functions (2in, 3in, 4in) ---------------------------


# st = '11100001'
# express = Boolean_dict[st]# ~d
# convertion = dict([('Or', 'And'),('And', 'Or')])
# func_converted = expr.expr(lf.multiple_replace(convertion, str(express[0])))
# bf = func_converted
# bf
#
# Opt_NOR_net = lf.permu_search_gens(bf,3,df3)
# Opt_NOR_net



# --------- generating 4inputs database -------
Num_total_fun, Boolean_dict = lf.Boolean_expresso_gens(4)
convertion = dict([('Or', 'And'),('And', 'Or')])
Database = []

for func, expression in Boolean_dict.items():
    # converted expression
    func_converted = expr.expr(lf.multiple_replace(convertion, str(expression[0])))
    Num_in = len(func_converted.inputs)
    # Final output of the NOR gates optimal network of the input boolean fucntion (1."And" "Or" Switched )
    Opt_NOR_net = lf.permu_search_gens(func_converted,4,df4)
    data = np.append(np.append([func,expression],Opt_NOR_net),Num_in)
    Database.append(list(data))

NOR4_database=np.array(Database)
NOR_df4 = pd.DataFrame({'Binary':NOR4_database[:,0],'Expresso_tts':np.array([i[0] for i in NOR4_database[:,1]]),'Num_gates':[i[0] for i in NOR4_database[:,2]],'Optimal_net':[i[1] for i in NOR4_database[:,2]], 'Permu_int':NOR4_database[:,3],'Num_in': NOR4_database[:,4]})
NOR_df4 = NOR_df4[NOR_df4.Num_in ==4]
NOR_df4.to_csv("NOR_df4.csv", index=False)



# --------- generating 3inputs database -------
Num_total_fun3, Boolean_dict3 = lf.Boolean_expresso_gens(3)
convertion = dict([('Or', 'And'),('And', 'Or')])
Database3 = []
for func, expression in Boolean_dict3.items():
    # converted expression
    func_converted = expr.expr(lf.multiple_replace(convertion, str(expression[0])))
    Num_in = len(func_converted.inputs)
    # Final output of the NOR gates optimal network of the input boolean fucntion (1."And" "Or" Switched )
    Opt_NOR_net = lf.permu_search_gens(func_converted,3,df3)
    data = np.append(np.append([func,expression],Opt_NOR_net),Num_in)
    Database3.append(list(data))


NOR3_database=np.array(Database3)
NOR_df3 = pd.DataFrame({'Binary':NOR3_database[:,0],'Expresso_tts':np.array([i[0] for i in NOR3_database[:,1]]),'Num_gates':[i[0] for i in NOR3_database[:,2]],'Optimal_net':[i[1] for i in NOR3_database[:,2]], 'Permu_int':NOR3_database[:,3],'Num_in': NOR3_database[:,4]})
NOR_df3 = NOR_df3[NOR_df3.Num_in ==3]
NOR_df3.to_csv("NOR_df3.csv", index=False)




# --------- generating 2inputs database -------
Num_total_fun2, Boolean_dict2 = lf.Boolean_expresso_gens(2)
convertion = dict([('Or', 'And'),('And', 'Or')])
Database2 = []
for func, expression in Boolean_dict2.items():
    # converted expression
    func_converted = expr.expr(lf.multiple_replace(convertion, str(expression[0])))
    Num_in = len(func_converted.inputs)
    # Final output of the NOR gates optimal network of the input boolean fucntion (1."And" "Or" Switched )
    Opt_NOR_net = lf.permu_search_gens(func_converted,2,df2)
    data = np.append(np.append([func,expression],Opt_NOR_net),Num_in)
    Database2.append(list(data))


NOR2_database=np.array(Database2)
NOR_df2 = pd.DataFrame({'Binary':NOR2_database[:,0],'Expresso_tts':np.array([i[0] for i in NOR2_database[:,1]]),'Num_gates':[i[0] for i in NOR2_database[:,2]],'Optimal_net':[i[1] for i in NOR2_database[:,2]], 'Permu_int':NOR2_database[:,3],'Num_in': NOR2_database[:,4]})
NOR_df2 = NOR_df2[NOR_df2.Num_in ==2]
NOR_df2.to_csv("NOR_df2.csv", index=False)



Full_database = [NOR_df2,NOR_df3,NOR_df4]
NOR_DATABASE = pd.concat(Full_database,ignore_index=True)
NOR_DATABASE.to_csv("NOR_DATABASE_filtered.csv", index=False)














# --------- ########## --------- #########-------------------   visualizations  --------- ########## --------- #########-------------------

#  ------ First reimport data ------
import pandas as pd
DB_import = pd.read_csv('NOR_DATABASE_filtered.csv')
# Add TN_functions column
from math import factorial
DB_import['TN_functions'] = DB_import['Num_in'].apply(lambda x : factorial(x))


# Import visualizations modules
import IPython
import altair as alt
def vegify(spec):
    IPython.display.display({
        'application/vnd.vegalite.v2+json': spec.to_dict()
    }, raw=True)
alt.data_transformers.enable('default', max_rows=None)
import seaborn as sns




DB_import
test = DB_import.iloc[0:400,:]

p2 = alt.Chart(DB_import.iloc[0:400,:]).mark_bar().encode(
    x='Num_gates',
    y='sum(TN_functions)',
    color=alt.Color('Num_in:N', scale=alt.Scale(scheme='set1'))
).configure(background="White").interactive()
vegify(p2)



p3 = alt.Chart(DB_import.iloc[0:400,:]).mark_bar().encode(
    x='Num_in',
    y='sum(Num_gates)',
    order=alt.Order('sum(Num_gates)', sort='ascending'),
    color=alt.Color('Num_gates:N', scale=alt.Scale(scheme='set1'))
).configure(background="White").interactive()
vegify(p3)



alt.Chart(DB_import.iloc[0:400,:]).mark_bar(opacity=0.7).encode(
    x='Num_gates:O',
    y=alt.Y('sum(TN_functions):Q', stack=None),
    color=alt.Color('Num_in:N', scale=alt.Scale(scheme='dark2')),
).configure(background="White")




# ---- normalized total number of functions with repect to each number of gate -----
# Relative abundance of different # of inputs for each min gates realization
relative_ab = alt.Chart(DB_import).mark_bar(opacity=0.7).encode(
    x=alt.X('Num_gates:O',title="# of min NOR gates"),
    y=alt.Y('sum(TN_functions):Q', stack="normalize", title="Fraction of total # of boolean functions"),
    color=alt.Color('Num_in:N', scale=alt.Scale(scheme='dark2'), title='Number of inputs'),
).properties(title = "Relative abundance of different # of inputs for each min gates realization"
            ).configure(background="White").interactive()

relative_ab.save('./Plots/relative_abundance.html')







alt.Chart(DB_import).mark_bar(opacity=0.7).encode(
    x='Num_in',
    y=alt.Y('sum(TN_functions):Q', stack="normalize"),
    # order=alt.Order('sum(Num_gates)', sort='ascending'),
    color=alt.Color('Num_in:N', scale=alt.Scale(scheme='set1')),
    column='Num_gates:N',
).properties(
    width=40
).configure(background="White").interactive()








# --- Plot3 :  Total number of functions histogram for each min NOR net in different # input cases (for p-class: just change Y->cout())
hist1 = alt.Chart(DB_import).mark_bar().encode(
    x=alt.X('Num_gates:O',title="# of min NOR gates"),
    y=alt.Y('sum(TN_functions):Q',title='Total # of boolean functions'),
    color=alt.Color('Num_gates:O',scale=alt.Scale(scheme='category20b')),
    column='Num_in:N'
).properties(title = "Total number of functions histogram for each min NOR net in different # input cases"
            ).configure(background="White").interactive()
hist1.save('./Plots/hist_gates_in.html')




# --- Plot4 :  Total number of functions histogram(discrete dot: easy to read) for each min NOR net in different # input cases (for p-class: just change Y->cout())
hist2 = alt.Chart(DB_import).mark_bar().encode(
    x=alt.X('Num_gates:O',title="# of min NOR gates"),
    y=alt.Y('sum(TN_functions):O',sort ='ascending',title='Total # of boolean functions'),
    color=alt.Color('Num_gates:O',scale=alt.Scale(scheme='category20b')),
    column='Num_in:N'
).properties(title = "Total number of functions square plot for each min NOR net in different # input cases"
            ).configure(background="White").interactive()
hist2.save('./Plots/hist_gates_in_square.html')












slider = alt.binding_range(min=1, max=14, step=1)
select_number_of_gates = alt.selection_single(name="NOR", fields=['Num_gates'], bind=slider)


slide_bar = alt.Chart(DB_import).mark_bar().encode(
    x=alt.X('Num_gates:N', axis=alt.Axis(title=None)),
    y=alt.Y('sum(TN_functions):O',title='Total # of boolean functions'),
    color=alt.Color('Num_in:N', scale=alt.Scale(scheme='set1')),
    column='Num_in:O'
).properties(title = "Total number of functions square plot for each min NOR net in different # input cases",
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
color=alt.Color('Num_gates:O', scale=alt.Scale(scheme='tableau20'), legend=alt.Legend(title="Number gates"))
# We create two selections:
# - a brush that is active on the top panel
# - a multi-click that is active on the bottom panel
brush = alt.selection_interval(encodings=['x'])
click = alt.selection_multi(encodings=['color'])


norm_bar = alt.Chart().mark_bar(opacity=0.7).encode(
    x=alt.X('Num_gates:O',title = 'minimun NOR gates'),
    y=alt.Y('sum(TN_functions):Q', stack="normalize", title = 'Total # of functions colored by inputs'),
    color=alt.Color('Num_in:O',legend=alt.Legend(title="Number inputs")),
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
    x=alt.X('Num_gates:N',title = 'minimun NOR gates'),
    y=alt.Y('sum(TN_functions):Q',title ='Total # of functions of all inputs'),
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
    y=alt.Y('sum(TN_functions):Q', stack="normalize")
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
f1=plt.figure(figsize=(12,10))
sns.violinplot(x="Num_in", y="Num_gates" ,bw= 0.1,scale="width",inner="quartile",data=DB_import)
# ax.xaxis.grid(True)
plt.xlabel("Number of inputs")
plt.ylabel("Number minimum NOR gates")
sns.despine(trim=True)
plt.tight_layout()
plt.axhline(y=7, color='r', linestyle=':')
plt.title("Distribution of min NOR gates(Categorical[1~14]) of functions with different inputs[2~4] ")
sns.set_context("poster")
plt.show()

f1.savefig("./Plots/violin_each_inputs.png", dpi = 600)




