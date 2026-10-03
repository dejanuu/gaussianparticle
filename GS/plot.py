import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

def read_obj(file_path):
    """
    Read an OBJ file and extract vertices and faces.
    
    Args:
        file_path: Path to the OBJ file
        
    Returns:
        vertices: numpy array of shape (n_vertices, 3)
        faces: list of lists containing vertex indices for each face
    """
    vertices = []
    faces = []
    
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#'):  # Skip comments
                continue
                
            values = line.split()
            if not values:
                continue
                
            if values[0] == 'v':  # Vertex
                vertices.append([float(values[1]), float(values[2]), float(values[3])])
                
            elif values[0] == 'f':  # Face
                # OBJ file indices start from 1, so we subtract 1
                # Also handle different face formats (v, v/vt, v/vt/vn)
                face = []
                for v in values[1:]:
                    vertex_index = int(v.split('/')[0]) - 1
                    face.append(vertex_index)
                faces.append(face)
    
    return np.array(vertices), faces

def plot_obj(file_path):
    """
    Read an OBJ file and plot it using matplotlib.
    
    Args:
        file_path: Path to the OBJ file
    """
    vertices, faces = read_obj(file_path)
    
    # Create a new figure and 3D axis
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Create mesh collection from the faces
    mesh = []
    for face in faces:
        # Get vertices for each face
        verts = [vertices[i] for i in face]
        mesh.append(verts)
    
    # Create a Poly3DCollection
    poly = Poly3DCollection(mesh, alpha=0.5, edgecolor='k', linewidth=0.5)
    poly.set_facecolor('lightblue')
    ax.add_collection3d(poly)
    
    # Auto-scale to the mesh size
    if len(vertices) > 0:
        x = vertices[:, 0]
        y = vertices[:, 1]
        z = vertices[:, 2]
        
        # Create equal aspect ratio
        max_range = np.array([x.max()-x.min(), y.max()-y.min(), z.max()-z.min()]).max() / 2.0
        mid_x = (x.max()+x.min()) / 2.0
        mid_y = (y.max()+y.min()) / 2.0
        mid_z = (z.max()+z.min()) / 2.0
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)
    
    # Add axis labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'3D Model: {file_path}')
    
    plt.tight_layout()
    plt.show()

# Call the function with your obj file
if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        obj_file = sys.argv[1]
    else:
        obj_file = "siris-output.obj"  # Default file or prompt for input
    
    plot_obj(obj_file)
