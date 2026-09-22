"""
Metabolomics comparison plots: noAra vs Ara, noAra vs iAAl, noAra vs Ara+iAAl, Ara vs Ara+iAAl
Produces, for each comparison: a heatmap (top variable compounds), a top-abundant bar plot,
and a top up/down (fold-change) bar plot.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style="whitegrid", font_scale=0.9)

DATA = "/home/claude/data"
OUT = "/home/claude/plots"

COMPARISONS = {
    "noAra_vs_Ara": {
        "file": f"{DATA}/EPI300_pCC1_ATF1_noAraAraOD.csv",
        "cols": ["EPI300_pCC1_ATF1_noAra", "EPI300_pCC1_ATF1_Ara"],
        "labels": ["noAra", "Ara"],
    },
    "noAra_vs_iAAl": {
        "file": f"{DATA}/EPI300_pCC1_ATF1_noAraiAAlOD.csv",
        "cols": ["EPI300_pCC1_ATF1_noAra", "EPI300_pCC1_ATF1_iAAl"],
        "labels": ["noAra", "iAAl"],
    },
    "noAra_vs_AraiAAl": {
        "file": f"{DATA}/EPI300_pCC1_ATF1_noAraAraiAAlOD.csv",
        "cols": ["EPI300_pCC1_ATF1_noAra", "EPI300_pCC1_ATF1_Ara_iAAl"],
        "labels": ["noAra", "Ara+iAAl"],
    },
    "Ara_vs_AraiAAl": {
        "file": f"{DATA}/EPI300_pCC1_ATF1_AraAraiAAlOD.csv",
        "cols": ["EPI300_pCC1_ATF1_Ara", "EPI300_pCC1_ATF1_AraAraiAAl"],
        "labels": ["Ara", "Ara+iAAl"],
    },
}

TOP_N_HEATMAP = 30
TOP_N_ABUNDANT = 20
TOP_N_UPDOWN = 15


def load_clean(entry):
    df = pd.read_csv(entry["file"])
    df = df.dropna(axis=1, how="all")  # drop stray empty trailing columns
    df = df.rename(columns={df.columns[0]: "Name"})
    df["Name"] = df["Name"].astype(str).str.strip()
    df = df[df["Name"].notna() & (df["Name"] != "nan") & (df["Name"] != "")]
    c1, c2 = entry["cols"]
    df[c1] = pd.to_numeric(df[c1], errors="coerce")
    df[c2] = pd.to_numeric(df[c2], errors="coerce")
    # collapse duplicate compound names by summing abundance (co-eluting hits)
    df = df.groupby("Name", as_index=False)[[c1, c2]].sum(min_count=1)
    # require detection in at least one condition
    df = df[(df[c1].notna()) | (df[c2].notna())]
    return df


def add_fc(df, c1, c2):
    # pseudocount = half the smallest nonzero, nonzero-NA value observed across both columns
    nonzero = pd.concat([df[c1], df[c2]])
    nonzero = nonzero[(nonzero > 0)]
    pseudo = nonzero.min() / 2 if len(nonzero) else 1.0
    a = df[c1].fillna(0) + pseudo
    b = df[c2].fillna(0) + pseudo
    df = df.copy()
    df["log2FC"] = np.log2(b / a)
    df["mean_abund"] = (df[c1].fillna(0) + df[c2].fillna(0)) / 2
    return df, pseudo


def plot_heatmap(df, entry, name):
    c1, c2 = entry["cols"]
    l1, l2 = entry["labels"]
    sub = df.reindex(df[[c1, c2]].fillna(0).sum(axis=1).sort_values(ascending=False).index)
    sub = sub.head(TOP_N_HEATMAP)
    mat = np.log10(sub[[c1, c2]].fillna(0) + 1)
    mat.index = sub["Name"].str.slice(0, 45)
    mat.columns = [l1, l2]

    fig, ax = plt.subplots(figsize=(6.5, max(6, 0.3 * len(mat))))
    sns.heatmap(mat, cmap="viridis", linewidths=0.4, linecolor="white",
                cbar_kws={"label": "log10(peak abundance + 1)"}, ax=ax,
                annot=False)
    ax.set_title(f"Top {len(mat)} most abundant compounds\n{l1} vs {l2}", fontsize=11)
    ax.set_ylabel("")
    plt.tight_layout()
    fig.savefig(f"{OUT}/{name}_heatmap.png", dpi=200)
    plt.close(fig)


def plot_top_abundant(df, entry, name):
    c1, c2 = entry["cols"]
    l1, l2 = entry["labels"]
    sub = df.copy()
    sub["mean_abund"] = (sub[c1].fillna(0) + sub[c2].fillna(0)) / 2
    sub = sub.sort_values("mean_abund", ascending=False).head(TOP_N_ABUNDANT)
    sub = sub.iloc[::-1]

    plot_df = sub.melt(id_vars="Name", value_vars=[c1, c2], var_name="Condition", value_name="Abundance")
    plot_df["Condition"] = plot_df["Condition"].map({c1: l1, c2: l2})
    plot_df["Name"] = pd.Categorical(plot_df["Name"], categories=sub["Name"], ordered=True)

    plot_df["Abundance"] = np.log10(plot_df["Abundance"].fillna(0) + 1)

    fig, ax = plt.subplots(figsize=(9, max(6, 0.35 * len(sub))))
    sns.barplot(data=plot_df, y="Name", x="Abundance", hue="Condition", ax=ax, palette="Set2")
    ax.set_title(f"Top {len(sub)} most abundant compounds (ranked by mean) — {l1} vs {l2}", fontsize=11)
    ax.set_xlabel("log10(peak abundance + 1)")
    ax.set_ylabel("")
    ax.legend(title="Condition", loc="lower right", framealpha=0.9)
    plt.tight_layout()
    fig.savefig(f"{OUT}/{name}_top_abundant.png", dpi=200)
    plt.close(fig)


def plot_top_updown(df, entry, name, pseudo):
    c1, c2 = entry["cols"]
    l1, l2 = entry["labels"]
    sub = df.dropna(subset=["log2FC"]).copy()
    up = sub.sort_values("log2FC", ascending=False).head(TOP_N_UPDOWN)
    down = sub.sort_values("log2FC", ascending=True).head(TOP_N_UPDOWN)
    combo = pd.concat([up, down]).drop_duplicates(subset="Name").sort_values("log2FC")
    combo["color"] = np.where(combo["log2FC"] > 0, "#d9534f", "#4b7bec")
    combo["label"] = combo["Name"].str.slice(0, 45)

    fig, ax = plt.subplots(figsize=(9, max(6, 0.3 * len(combo))))
    ax.barh(combo["label"], combo["log2FC"], color=combo["color"])
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_title(f"Top up/down compounds: {l2} vs {l1}\n(log2 fold change, pseudocount={pseudo:.1f})", fontsize=11)
    ax.set_xlabel(f"log2FC ({l2} / {l1})")
    plt.tight_layout()
    fig.savefig(f"{OUT}/{name}_top_updown.png", dpi=200)
    plt.close(fig)


summary_rows = []
for name, entry in COMPARISONS.items():
    df = load_clean(entry)
    df, pseudo = add_fc(df, *entry["cols"])
    plot_heatmap(df, entry, name)
    plot_top_abundant(df, entry, name)
    plot_top_updown(df, entry, name, pseudo)
    summary_rows.append({"comparison": name, "n_compounds": len(df), "pseudocount": pseudo})
    print(f"{name}: {len(df)} compounds, pseudocount={pseudo:.3f}")

print(pd.DataFrame(summary_rows))
