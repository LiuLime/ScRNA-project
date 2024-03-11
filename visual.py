"""Tools for Visualization


2024/3/8 @Liu
"""
from sklearn.datasets import load_digits
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns

data=load_digits()