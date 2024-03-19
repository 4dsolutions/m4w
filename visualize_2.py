"""
Another AI endeavor
Copied from:  https://github.com/docxology/ivm-xyz/

The code below bears a family resemblance to code in my repos, including this one.
However AI is generating code from other sources and in this example is using 
threads, which I haven't used on my end for this kind of image file generating.

The GIFs this code produces are pretty huge, like 13MB per file, each one several
rotations of a polyhedron as depicted using matplotlib's 3d toolkit. Once again,
this is not an approach I've taken, even though I'm familiar with matplotlib. My
graphical output has been via raytracing with povray, infrequently in Blender, or
using VPython for realtime animation. I cut my teeth on VRML back in the day, and
POV-Ray has always been a favorite.

The AI-generated dodecahedron is a mess in this version. The others look OK.
No icosahedron was included.

Where to find more information about these AI experiments is in this Coda:
https://coda.io/d/Math4Wisdom_d0SvdI3KSto/ivm-xyz_suqdu#_luR7B

Warning: the code below takes minutes on one CPU and creates gigantic GIFs.
I'm not saying you shouldn't run it, just give it time to finish, don't 
conclude it's hung, as I did the first time I ran it.
"""

import matplotlib
matplotlib.use('Agg')  # Use a non-interactive backend to avoid RuntimeError with tkinter
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import os
import sys
import imageio
from concurrent.futures import ThreadPoolExecutor

# User-specific hardcoded maximum number of threads
MAX_THREADS = 8

class PolyhedronPlotter:
    def __init__(self, output_folder="images"):
        self.output_folder = output_folder
        os.makedirs(self.output_folder, exist_ok=True)
        self.polyhedrons = []  # Initialize an empty list to store polyhedron dictionaries

    def add_polyhedron(self, vertices, faces, title, file_name):
        """
        Adds a polyhedron to the list of polyhedrons to be animated.
        
        :param vertices: A list of tuples, each representing the x, y, z coordinates of a vertex.
        :param faces: A list of tuples, each representing indices of vertices that form a face.
        :param title: Title of the polyhedron.
        :param file_name: Name of the file to save the animation as.
        """
        polyhedron = {
            "vertices": vertices,
            "faces": faces,
            "title": title,
            "file_name": file_name
        }
        self.polyhedrons.append(polyhedron)

    def animate_polyhedron(self, save=False):
        """
        Creates an animated GIF of each polyhedron in the list by improving aesthetics and adding vertex labels.
        Optionally saves the animation as a GIF file. Now parallelized across as many CPU threads as possible, 
        up to a user-specific hardcoded maximum.

        :param save: Boolean indicating whether to save the plot as a GIF file.
        """
        num_threads = min(len(self.polyhedrons), MAX_THREADS)

        def plot_and_save(polyhedron):
            vertices = polyhedron["vertices"]
            faces = polyhedron["faces"]
            title = polyhedron["title"]
            file_name = polyhedron["file_name"]

            images = []
            for angle in range(0, 360, 2):
                fig = plt.figure(figsize=(10, 8))
                ax = fig.add_subplot(111, projection='3d')

                vtx = [[vertices[i] for i in face] for face in faces]

                poly = Poly3DCollection(vtx, facecolors='skyblue', linewidths=0.5, edgecolors='darkblue', alpha=0.5)
                ax.add_collection3d(poly)

                for i, (x, y, z) in enumerate(vertices):
                    ax.scatter(x, y, z, color="darkred", s=100, edgecolors='black', zorder=5)
                    ax.text(x, y, z, f'V{i+1} ({x}, {y}, {z})', color='black')

                ax.set_xlabel('X Axis', fontsize=12)
                ax.set_ylabel('Y Axis', fontsize=12)
                ax.set_zlabel('Z Axis', fontsize=12)

                ax.set_xlim([min(v[0] for v in vertices)-1, max(v[0] for v in vertices)+1])
                ax.set_ylim([min(v[1] for v in vertices)-1, max(v[1] for v in vertices)+1])
                ax.set_zlim([min(v[2] for v in vertices)-1, max(v[2] for v in vertices)+1])

                ax.view_init(30, angle)
                plt.title(title, fontsize=14, fontweight='bold')
                plt.tight_layout()

                # Convert the plot to an image array
                fig.canvas.draw()
                image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
                image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))

                images.append(image)

                plt.close(fig)

            imageio.mimsave(os.path.join(self.output_folder, file_name), images, fps=20)

        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            executor.map(plot_and_save, self.polyhedrons)

if __name__ == "__main__":
    plotter = PolyhedronPlotter()
    # Example usage
    platonic_solids = [
        {
            "vertices": [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)],  # Tetrahedron
            "faces": [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)],
            "title": "Tetrahedron",
            "file_name": "tetrahedron.gif"
        },
        {
            "vertices": [(0, 0, 1), (1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, -1, 0)],  # Pyramid (Square Base)
            "faces": [(0, 1, 2), (0, 2, 3), (0, 3, 4), (0, 4, 1), (1, 2, 3, 4)],
            "title": "Pyramid",
            "file_name": "pyramid.gif"
        },
        {
            "vertices": [(0, 0, -1), (0, 0, 1), (1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, -1, 0)],  # Octahedron
            "faces": [(0, 2, 3), (0, 3, 4), (0, 4, 5), (0, 5, 2), (1, 2, 3), (1, 3, 4), (1, 4, 5), (1, 5, 2)],
            "title": "Octahedron",
            "file_name": "octahedron.gif"
        },
        {
            "vertices": [(-1, -1, -1), (-1, -1, 1), (-1, 1, -1), (-1, 1, 1), (1, -1, -1), (1, -1, 1), (1, 1, -1), (1, 1, 1)],  # Cube
            "faces": [(0, 1, 3, 2), (4, 6, 7, 5), (0, 2, 6, 4), (1, 5, 7, 3), (0, 4, 5, 1), (2, 3, 7, 6)],
            "title": "Cube",
            "file_name": "cube.gif"
        },
        {
            "vertices": [(0, 0, 1.176), (1.051, 0, 0.526), (0.324, 1, 0.525), (-0.851, 0.618, 0.526), (-0.851, -0.618, 0.526), (0.325, -1, 0.526), (0.851, 0.618, -0.526), (0.851, -0.618, -0.526), (-0.325, 1, -0.526), (-1.051, 0, -0.526), (-0.325, -1, -0.526), (0, 0, -1.176)],  # Dodecahedron
            "faces": [(0, 1, 2, 3, 4), (0, 4, 5, 6, 1), (1, 6, 7, 8, 2), (2, 8, 9, 3), (3, 9, 10, 4), (5, 4, 10, 11, 6), (6, 11, 7), (7, 11, 10, 9, 8), (0, 2, 1), (5, 4, 0)],
            "title": "Dodecahedron",
            "file_name": "dodecahedron.gif"
        }
    ]
    for solid in platonic_solids:
        plotter.add_polyhedron(**solid)
    # Add more polyhedrons as needed
    plotter.animate_polyhedron(save=True)