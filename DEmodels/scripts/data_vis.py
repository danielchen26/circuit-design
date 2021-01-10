import pandas as pd
import altair as alt
alt.data_transformers.disable_max_rows()
# pd.options.display.html.table_schema = True
df2 = pd.read_csv("2Bits_DB.csv")
df2
alt.Chart(df2).mark_point().encode(
    x='K',
    y='n',
    # color='δ',
)
