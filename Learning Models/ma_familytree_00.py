# -*- coding: utf-8 -*-
"""
Filename: ma_familytree_00.py
Date created: 2022/11/09
@author: Mehvish Ashiq
Purpose: 
Steps: 
https://www.delftstack.com/howto/python-pandas/pandas-family-tree/
Maybe: https://stackoverflow.com/questions/66823677/how-to-draw-a-family-tree-from-a-pandas-dataframe
"""
import pandas as pd
import numpy as np
from graphviz import Digraph

rawdf = pd.read_csv("./data.csv", keep_default_na=False)
element1 = rawdf[["id", "mid"]]
element2 = rawdf[["id", "fid"]]
element1.columns = ["Child", "ParentID"]
element2.columns = element1.columns

element = pd.concat([element1, element2])
element.replace("", np.nan, regex=True, inplace=True)
t = pd.DataFrame({"tmp": ["no_entry" + str(i) for i in range(element.shape[0])]})
element["ParentID"].fillna(t["tmp"], inplace=True)
df = element.merge(rawdf, left_index=True, right_index=True, how="left")

df["name"] = df[df.columns[4:6]].apply(
    lambda x: " ".join(x.dropna().astype(str)), axis=1
)
df = df.drop(["Child", "fid", "mid", "first_name", "last_name"], axis=1)
df = df[["id", "name", "gender", "dob", "dod", "birth_place", "job", "ParentID"]]

f = Digraph(
    "neato",
    format="pdf",
    encoding="utf8",
    filename="data",
    node_attr={"color": "lightblue2", "style": "filled"},
)
f.attr("node", shape="box")

for index, record in df.iterrows():
    f.edge(str(record["ParentID"]), str(record["id"]), label="")
f.view()

f = Digraph(
    "neato",
    format="jpg",
    encoding="utf8",
    filename="detailed_data",
    node_attr={"style": "filled"},
    graph_attr={"concentrate": "true", "splines": "ortho"},
)
f.attr("node", shape="box")

for index, row in df.iterrows():
    f.node(
        row["id"],
        label=row["name"]
        + "\n"
        + row["job"]
        + "\n"
        + str(row["dob"])
        + "\n"
        + row["birth_place"]
        + "\n"
        + str(row["dod"]),
        _attributes={
            "color": "lightpink"
            if row["gender"] == "F"
            else "lightblue"
            if row["gender"] == "M"
            else "lightgray"
        },
    )

for index, row in df.iterrows():
    f.edge(str(row["ParentID"]), str(row["id"]), label="")
f.view()