# -*- coding: utf-8 -*-
"""
Filename: ma_familytree.py
Date created: 2024/04/17, Wed, 20:50:37 (UTC+8)
@author: Mehvish Ashiq, Lio Hong
Purpose: Use pandas library and graphviz to create a family tree.
Steps: 
https://www.delftstack.com/howto/python-pandas/pandas-family-tree/
Maybe: https://stackoverflow.com/questions/66823677/how-to-draw-a-family-tree-from-a-pandas-dataframe
data.csv copied from table provided in first link.
"""
import pandas as pd
import numpy as np
from graphviz import Digraph

# To adjust the dataframe appearance
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 20)
pd.set_option("display.width", 200)
pd.set_option('display.expand_frame_repr', False)

rawdf = pd.read_csv("./data.csv", keep_default_na=False)
# Separate by each parent, which is also the edge list format.
element1 = rawdf[["id", "mid"]]
element2 = rawdf[["id", "fid"]]
element1.columns = ["Child", "ParentID"]
element2.columns = element1.columns

# To bottom of element1, add rows of element2.
element = pd.concat([element1, element2])
# Prepare df to use df.fillna().
element.replace("", np.nan, regex=True, inplace=True)
# Creates a 1-col df of no_entry+i strings.
t = pd.DataFrame({"tmp": ["no_entry" + str(i) for i in range(element.shape[0])]})
element["ParentID"].fillna(t["tmp"], inplace=True)
df = element.merge(rawdf, left_index=True, right_index=True, how="left")

df["name"] = df[df.columns[4:6]].apply(
    lambda x: " ".join(x.dropna().astype(str)), axis=1
)
df = df.drop(["Child", "fid", "mid", "first_name", "last_name"], axis=1)
df = df[["id", "name", "gender", "dob", "dod", "birth_place", "job", "ParentID"]]

def basic():
    f = Digraph(
        "neato",
        format="pdf",
        encoding="utf8",
        filename="data",
        node_attr={"color": "lightblue2", "style": "filled"},
    )
    f.attr("node", shape="box")
    
    for index, record in df.iterrows():
        # Source to destination.
        f.edge(str(record["ParentID"]), str(record["id"]), label="")
    f.view()


def bbasic():
    return None


def tinted():
    f = Digraph(
        "neato",
        format="jpg",
        encoding="utf8",
        filename="detailed_data",
        node_attr={"style": "filled"},
        graph_attr={"concentrate": "true", "splines": "ortho"},
    )
    f.attr("node", shape="box")
    
    # Index not used directly, but needed for syntax.
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
                # This color never appears, but needed for syntax.
                else "lightgray",
                "style": "filled"
            },
        )
    
    for index, row in df.iterrows():
        f.edge(str(row["ParentID"]), str(row["id"]), label="")
    f.view()
    

# This might be foolhardy but I wonder if I can handle individuals without parents.
def ttinted():
    f = Digraph(
        "neato",
        format="jpg",
        encoding="utf8",
        filename="clean_data",
        node_attr={"style": "invis"},
        graph_attr={"concentrate": "true", "splines": "ortho"},
        edge_attr={"style": "invis"}
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
                else "lightgray",
                "style": "filled"
            },
        )
    
    '''
    Removing this breaks the tree in 3 ways:
    1. Flattens entire tree.
    2. Removes all arrows.
    3. Removes all parents not present as children.
    '''
    for index, row in df.iterrows():
        f.edge(str(row["ParentID"]), str(row["id"]), label="",
            _attributes={
               "style": "filled"
               if row["ParentID"] in df.id.values
               else "invis"})
    f.view()
    return f
f = ttinted()