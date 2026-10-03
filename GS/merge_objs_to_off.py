import os

def read_obj(file_path):
    vertices = []
    faces = []

    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            values = line.strip().split()
            if not values:
                continue
            if values[0] == 'v':
                vertices.append(list(map(float, values[1:4])))
            elif values[0] == 'f':
                # Only handle 'f v1 v2 v3' (no texture/normal)
                face = [int(part.split('/')[0]) - 1 for part in values[1:]]
                faces.append(face)

    return vertices, faces

def write_off(file_path, vertices, faces):
    with open(file_path, 'w') as f:
        f.write("OFF\n")
        f.write(f"{len(vertices)} {len(faces)} 0\n")
        for v in vertices:
            f.write(f"{v[0]} {v[1]} {v[2]}\n")
        for face in faces:
            f.write(f"{len(face)} {' '.join(map(str, face))}\n")
    print(f"✅ Merged OFF file written to: {file_path}")

def merge_obj_to_off(obj_files, output_off_path):
    all_vertices = []
    all_faces = []
    vertex_offset = 0

    for file in obj_files:
        vertices, faces = read_obj(file)
        all_vertices.extend(vertices)

        # Offset face indices based on current total vertex count
        offset_faces = [[idx + vertex_offset for idx in face] for face in faces]
        all_faces.extend(offset_faces)

        vertex_offset += len(vertices)

    write_off(output_off_path, all_vertices, all_faces)

# Example usage
if __name__ == "__main__":
    import sys
    if len(sys.argv) > 2:
        obj_files = sys.argv[1:-1]  # all except last
        output_off = sys.argv[-1]   # last argument
    else:
        print("Usage: python merge_objs_to_off.py model1.obj model2.obj ... output.off")
        sys.exit(1)

    merge_obj_to_off(obj_files, output_off)
