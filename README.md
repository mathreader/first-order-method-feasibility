# Subgradient Method on Linear Programming Feasibilty Problems

This repository contains code for the report [First-order Methods on Large-scale Convex Optimization Problems](https://mathreader.github.io/files/report_first_order_method.pdf)

The main testing file to run is test_subgrad.m. This contains comparison between the following methods:
1. Subgradient method with 1/n step size.
1. Subgradient method with fixed step size eps/|g_k|^2.
1. Subgradient method with Polyak's step size.
1. Subgradient method with Polyak first, then fixed step size.
1. Subgradient method with restart.
1. Perceptron method with fixed step size eps/|g_k|^2.
1. Perceptron method with Polyak's step size.

The other testing files focus on analyzing the effect of different parameters on the restart method. The results for those experiments are contained in the ./data/ folder, and can be accessed by running data_manipulation.m within the folder.

