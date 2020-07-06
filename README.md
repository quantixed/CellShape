# CellShape

[![DOI](https://zenodo.org/badge/202709292.svg)](https://zenodo.org/badge/latestdoi/202709292)

Analysis of cell shape from 2D outlines in IGOR Pro.

[**Analysis**](#analysis) | [**Outputs**](#outputs) | [**Input Data**](#input-data)

## Analysis

A number of cell shape descriptors are analyzed.
Note that the intended use of `CellShape` for analysis of cell outlines from single images of different cells grouped by experimental condition, i.e. it is not (currently) intended for use with movie-type data.

Briefly, **violin plots** where dots represent individual cells and the marker indicates the mean, show:

- minor axis length
- major axis length
- cell perimeter
- cell area (mask area is an equivalent metric)
- aspect ratio, minor:major axes
- symmetry, ratio of the cell area to the area of a the cell footprint reflected on its major axis
- convex area, area of a convex hull enclosing the cell perimeter
- solidity, ratio of cell area to convex area
- extent, ratio of cell area to the area of a bounding box
- maximum distance to the perimeter from the cell center
- minimum distance to the perimeter from the cell center

A **plot** of minor axis length vs major axis length.
Each cell in the dataset is represented by a bubble, the size of which corresponds to that cell's area.

An **image quilt** to show the dataset.

## Outputs

*Example data is available in this repo together with example outputs.*

A series of plots that compare the experimental conditions are generated.

![img](img/output1.png?raw=true "image")

Additionally, an image quilt of cell outlines is generated, this is maximized to show the most data possible while making same-sized quilts.

![img](img/output2.png?raw=true "image")

## Input data

Cell outlines are hand drawn by a user blind to the conditions of the experiment.
We construct IMOD models - 1 model per image, multiple contours (outlines) per model.
For input into Igor we need to convert each model to a text file using `model2point`.
For example:

```bash
model2point -fl -ob -z example_Model_IMOD example.txt
```

In the shell, a directory of models can be converted using:

```bash
find . -type f ! -name "*.*" |
 while IFS= read file_name; do
 model2point -fl -ob -z "$file_name" "${file_name##*\/}.txt"
 done
```

If using another package to generate the outlines, the text files need to be space delimited, with no header.
Showing `model contour x y z`.
For example:

```
     0     0      104.00     1229.00        0.00
     0     0      100.00     1225.00        0.00
     0     0      100.00     1219.00        0.00
```

## Running the code

1. Save `CellShape.ipf` in *Wavemetrics/Igor Pro 8 User Files/User Procedures*.
2. Open in Igor and compile.
3. Run using Macros>Cell Shape...
4. Select the folder of point files to be analysed.

Igor will determine the conditions, as long as the files are logically named (see ipf for details).
Igor will ask in which order you'd like the outputs organised.

![img](img/inOrder.png?raw=true "image")

Next, you have the opportunity to rename conditions for the plots to make sense.

![img](img/aliasCheck.png?raw=true "image")

Click **Do It** and Igor will do the rest!