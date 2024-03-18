"""
This code was developed using AI tools as an experiment that started with 
the 4dsolutions repos on Github. 

More details here:  https://coda.io/d/Math4Wisdom_d0SvdI3KSto/ivm-xyz_suqdu#_luR7B

"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import os
import sys

class PolyhedronPlotter:
    def __init__(self, output_folder="images"):
        self.output_folder = output_folder

    def plot_polyhedron(self, vertices, faces, title="Polyhedron Visualization", save=False, file_name="polyhedron.png"):
        """
        Enhances the visualization of a polyhedron given its vertices and faces by improving aesthetics and adding vertex labels.
        Optionally saves the plot as a PNG file.

        :param vertices: A list of tuples, each representing the x, y, z coordinates of a vertex.
        :param faces: A list of tuples, each representing indices of vertices that form a face.
        :param title: Title of the plot.
        :param save: Boolean indicating whether to save the plot as a PNG file.
        :param file_name: Name of the file to save the plot as. Defaults to 'polyhedron.png'.
        """
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')

        vtx = [[vertices[i] for i in face] for face in faces]

        poly = Poly3DCollection(vtx, facecolors='skyblue', linewidths=0.5, edgecolors='darkblue', alpha=0.5)
        ax.add_collection3d(poly)

        for i, (x, y, z) in enumerate(vertices):
            ax.scatter(x, y, z, color="darkred", s=100, edgecolors='black', zorder=5)
            ax.text(x, y, z, f'  V{i+1} ({x}, {y}, {z})', color='black')

        ax.set_xlabel('X Axis', fontsize=12)
        ax.set_ylabel('Y Axis', fontsize=12)
        ax.set_zlabel('Z Axis', fontsize=12)

        ax.set_xlim([min(v[0] for v in vertices)-1, max(v[0] for v in vertices)+1])
        ax.set_ylim([min(v[1] for v in vertices)-1, max(v[1] for v in vertices)+1])
        ax.set_zlim([min(v[2] for v in vertices)-1, max(v[2] for v in vertices)+1])

        plt.title(title, fontsize=14, fontweight='bold')
        plt.tight_layout()

        if save:
            if not os.path.exists(self.output_folder):
                os.makedirs(self.output_folder)
            plt.savefig(f"{self.output_folder}/{file_name}")
        else:
            plt.show()
        plt.close(fig)

if __name__ == "__main__":
    plotter = PolyhedronPlotter()
    polyhedron_list = [
        {
            "vertices": [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)],
            "faces": [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)],
            "title": "Tetrahedron 1",
            "file_name": "tetrahedron_1.png"
        },
        {
            "vertices": [(0, 0, 0), (2, 0, 0), (0, 2, 0), (0, 0, 2)],
            "faces": [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)],
            "title": "Tetrahedron 2",
            "file_name": "tetrahedron_2.png"
        }
    ]
    for polyhedron in polyhedron_list:
        plotter.plot_polyhedron(polyhedron["vertices"], polyhedron["faces"], title=polyhedron["title"], save=True, file_name=polyhedron["file_name"])