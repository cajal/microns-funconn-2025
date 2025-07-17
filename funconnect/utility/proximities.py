from pathlib import Path
import pandas as pd
import numpy as np
from scipy.spatial import KDTree
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from typing import Optional, Tuple, Sequence, Dict, Any

nm_per_vx = np.array([4, 4, 40]) 
seg_bounds_vx = np.array([[52770, 60616, 14850], [437618, 322718, 27858]]) # from graphene://https://minnie.microns-daf.com/segmentation/table/minnie3_v1
seg_bounds_nm = seg_bounds_vx * nm_per_vx

def load_sample_skeleton(
    nucleus_id,
    compartment,
    src_dir=None
):
    if src_dir is None:
        src_dir = Path(__file__).parent.parent.parent / 'methods' / 'sample_data' / 'skeletons'
    else:
        src_dir = Path(src_dir)
    lookup_df = pd.read_csv(src_dir / 'skeleton_lookup.csv')
    filename = lookup_df.loc[
        (lookup_df['nucleus_id'] == nucleus_id) &
        (lookup_df['compartment'] == compartment),
        'filename'
    ].values[0]
    skeleton_path = src_dir / filename
    return np.load(skeleton_path)


def convert_skeleton_to_nodes_edges(
    skeleton,
    verbose=False,
    ):
    """
    Convert a skeleton representation into unique nodes and edges.

    This function takes a skeleton represented as an array of 3D line segments
    (pairs of vertices) and returns a deduplicated set of nodes and an array
    of edges defined by indices into this node list.

    Parameters
    ----------
    skeleton : np.ndarray of shape (N, 2, 3)
        An array representing N line segments in 3D space. Each segment is
        defined by two 3D points.

    verbose : bool, optional
        If True, prints debugging information during processing. Default is False.

    Returns
    -------
    unique_rows : np.ndarray of shape (M, 3)
        Array of unique 3D coordinates representing the deduplicated set of nodes
        in the skeleton.

    reshaped_indices : np.ndarray of shape (N, 2)
        Array of index pairs into `unique_rows`, representing the edges between nodes.
    """
    
    all_skeleton_vertices = skeleton.reshape(-1,3)
    unique_rows,indices = np.unique(all_skeleton_vertices,return_inverse=True,axis=0)

    #need to merge unique indices so if within a certain range of each other then merge them together
    reshaped_indices = indices.reshape(-1,2)
    
    return unique_rows,reshaped_indices


def discretize_skeleton(skeleton, max_length, return_mapping=False):
    """
    Discretize the provided skeleton based on the specified maximum length.
    
    Parameters:
    - skeleton (ndarray): Array of shape (N, 2, 3) representing N edges with start and end points.
    - max_length (float): The maximum allowed length for any segment in the discretized skeleton.
    - return_mapping (bool): If True, a mapping array will be returned that would allow the original skeleton to be recovered. Default is False.
    
    Returns:
    - ndarray: Discretized skeleton based on the given parameters.
    """
    start_points, end_points = skeleton[:, 0], skeleton[:, 1]
    segment_diffs = end_points - start_points
    segment_lengths = np.linalg.norm(segment_diffs, axis=1)
    
    num_segments = np.ceil(segment_lengths / max_length).astype(int)
    num_segments[segment_lengths <= max_length] = 1

    if np.all(num_segments == 1):
        if return_mapping:
            return skeleton, np.arange(len(skeleton))
        return skeleton
    
    adjusted_diffs = segment_diffs / num_segments[:, None]
    segmented_diffs = np.repeat(adjusted_diffs, num_segments, axis=0)

    repeated_starts = np.repeat(start_points, num_segments, axis=0)
    increments = np.hstack([np.arange(n) for n in num_segments])
    interpolated_points = repeated_starts + increments[:, None] * segmented_diffs
    
    total_points = 0
    discretized_segments = []
    for idx, splits in enumerate(num_segments):
        segment_points = interpolated_points[total_points:total_points + splits]
        
        if splits == 1:
            segment_points = skeleton[idx]
        else:
            segment_points = np.vstack([segment_points, end_points[idx]])
        
        discretized_segments.append(segment_points)
        total_points += splits

    mapping = []
    for idx, splits in enumerate(num_segments):
        mapping.extend([idx] * splits)
    mapping = np.array(mapping)

    discretized_skeleton = np.vstack([np.array((segment[:-1], segment[1:])).transpose(1, 0, 2) for segment in discretized_segments])

    if return_mapping:
        return discretized_skeleton, mapping
    return discretized_skeleton


def discretized_to_original_skeleton(discretized_skeleton, mapping):
    """
    Convert a discretized skeleton back to its original form using the provided mapping.
    
    Parameters:
    - discretized_skeleton (ndarray): Array of shape (M, 2, 3) representing edges of a discretized skeleton.
    - mapping (ndarray): Mapping array of shape (M,) representing original edge index for each discretized edge.
    
    Returns:
    - ndarray: Original skeleton edges.
    """
    
    _, segment_counts = np.unique(mapping, return_counts=True)
    
    original_edges = []
    idx = 0
    for segment_count in segment_counts:
        start_point = discretized_skeleton[idx][0]
        end_point = discretized_skeleton[idx + segment_count - 1][1]
        idx += segment_count
        original_edges.append((start_point, end_point))
    
    return np.stack(original_edges)


def filter_edges(edges, vertices_inds_subset):
    """
    Find and return the subset of edges where at least one node is present in vertices_inds_subset.

    :param edges: List of lists or array of shape (N, 2) representing edges by their vertex indices.
    :param vertices_inds_subset: List or array of shape (M,) containing vertex indices to filter edges against.
    :return: Array of shape (P, 2) where P ≤ N, representing the filtered edges.
    """
    
    inds_set = set(vertices_inds_subset)
    filtered_edges = [edge for edge in edges if edge[0] in inds_set or edge[1] in inds_set]
    if not filtered_edges:
        return np.array([], dtype=int).reshape(0, 2)
    else:
        return np.stack(filtered_edges)


def compute_skeletal_length(vertices, edges, vertices_inds_subset=None):
    """
    Compute the total length of the skeletal structure given vertices and edges. 
    If a subset of vertex indices is provided, it will calculate the length only for 
    the edges related to that subset.

    :param vertices: N x 3 array-like object containing the coordinates of each vertex.
    :param edges: N x 2 array-like object of edge pairs, where each edge is represented by the indices of its endpoints.
    :param vertices_inds_subset: (Optional) A list of vertex indices for which to filter edges with. 
                                 If not provided, the function computes the skeletal length for all edges.
    :return: The total length of the skeletal structure (or the subset if vertices_inds_subset is provided).
    """
    vertices = np.array(vertices)
    edges = np.array(edges)
    
    if vertices_inds_subset is not None:
        edges = filter_edges(edges, vertices_inds_subset)

    assert vertices.shape[1] == 3, f"Invalid vertices shape: {vertices.shape}"
    assert edges.shape[1] == 2, f"Invalid edge shape: {edges.shape}"
    
    skeleton = vertices[edges]

    # Calculate the differences for each edge's endpoints
    diffs = skeleton[:, 1, :] - skeleton[:, 0, :]
    
    # Compute the norm of each difference (edge length)
    edge_lengths = np.linalg.norm(diffs, axis=1)
    
    return edge_lengths.sum()


def compute_proximities(verts1, verts2, radius):
    """
    Compute pairwise proximities between two sets of 3D vertices using KDTree.

    For each point in `verts1`, this function finds all points in `verts2` that 
    lie within a specified Euclidean distance (`radius`). The function returns 
    the matched vertex pairs, their original indices, and the corresponding distances.

    Parameters
    ----------
    verts1 : np.ndarray of shape (N, 3)
        First set of 3D points (typically a surface mesh or point cloud).
    
    verts2 : np.ndarray of shape (M, 3)
        Second set of 3D points to compare with `verts1`.

    radius : float
        Distance threshold for proximity. Only pairs of points within this
        distance are considered proximal.

    Returns
    -------
    dict
        A dictionary containing:
        - 'verts1_prx': np.ndarray of shape (K, 3)
            Subset of points from `verts1` that are within `radius` of at least one point in `verts2`.
        - 'verts2_prx': np.ndarray of shape (K, 3)
            Corresponding points in `verts2` matched to `verts1_prx`.
        - 'verts1_inds_prx': np.ndarray of shape (K,)
            Indices of `verts1_prx` in the original `verts1` array.
        - 'verts2_inds_prx': np.ndarray of shape (K,)
            Indices of `verts2_prx` in the original `verts2` array.
        - 'dists_prx': np.ndarray of shape (K,)
            Rounded Euclidean distances between each matched pair of points.

        If no proximal pairs are found, an empty dictionary is returned.
    """

    # initialize result variables
    verts1_prx = []
    verts2_prx = []
    verts1_inds_prx = []
    verts2_inds_prx = []
    dists_prx = []

    # compute proximities
    tree1 = KDTree(verts1)
    tree2 = KDTree(verts2)
    indices2 = tree1.query_ball_tree(tree2, r=radius)

    # extract proximities
    for ind1, inds2_group in enumerate(indices2):
        if len(inds2_group) > 0:
            for ind2 in inds2_group:
                verts1_inds_prx.append(ind1)
                verts2_inds_prx.append(ind2)

    if len(verts1_inds_prx) == 0:
        return {}

    verts1_inds_prx = np.array(verts1_inds_prx)
    verts2_inds_prx = np.array(verts2_inds_prx)
    verts1_prx = verts1[verts1_inds_prx]
    verts2_prx = verts2[verts2_inds_prx]
    dists_prx = np.linalg.norm(verts1_prx - verts2_prx, axis=1).round().astype(int)

    return {
            'verts1_prx': verts1_prx,
            'verts2_prx': verts2_prx,
            'verts1_inds_prx': verts1_inds_prx,
            'verts2_inds_prx': verts2_inds_prx,
            'dists_prx': dists_prx,
    }


### PLOTTING FUNCTIONS 

def add_skeleton(ax, vertices, edges, xdim=0, ydim=1, edge_skip=1, *args, **kwargs):
    """
    Adds a skeleton to the current plot.

    Parameters:
    - vertices: Nx3 array of vertex coordinates.
    - edges: Sequence of (i, j) index pairs into vertices array.
    - xdim: Dimension index for x-coordinates (default: 0).
    - ydim: Dimension index for y-coordinates (default: 1).
    - edge_skip: Plot every nth edge (default: 1, plots all edges).
    - *args, **kwargs: Additional arguments passed to the plot function.
    Returns:
    - line: The matplotlib Line2D object representing the skeleton.
    """

    for edge in edges[::edge_skip]:
        line, = ax.plot(*vertices[edge].T, *args, **kwargs)
    return line


def add_scale_bar(ax, length='auto', fraction=0.1, location='lower left', linewidth=2, text_offset=5, unit_label='$\mu m$', **kwargs):
    """
    Adds a horizontal scale bar to the given axes.

    Parameters:
    - ax: matplotlib axes object
    - length: float or 'auto'. If 'auto', uses fraction * x-range.
    - fraction: used only if length='auto'; determines the scale bar size.
    - location: position of the scale bar ('lower left' supported)
    - linewidth: width of the scale bar line
    - text_offset: vertical offset in points between scale bar and text label
    - kwargs: extra arguments passed to Line2D
    """
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    x_range = xlim[1] - xlim[0]
    y_range = ylim[1] - ylim[0]

    # Determine bar length
    if length == 'auto':
        bar_length = x_range * fraction
    else:
        bar_length = float(length)

    # Location of the scale bar
    if location == 'lower left':
        x_start = xlim[0] + 0.1 * x_range
        y_start = ylim[0] + 0.1 * y_range
    else:
        raise NotImplementedError("Only 'lower left' location is implemented.")

    # Add scale bar line
    bar = mlines.Line2D([x_start, x_start + bar_length], [y_start, y_start],
                        color='black', linewidth=linewidth, **kwargs)
    ax.add_line(bar)

    # Convert text_offset from points to data coordinates
    text_offset_data = text_offset * ax.figure.dpi / 72 / ax.get_window_extent().height * y_range

    # Add text label centered below the bar
    ax.text(x_start + bar_length / 2, y_start - text_offset_data,
            f'{bar_length:.0f} {unit_label}', ha='center', va='top')
    

def expand_bbox_to_aspect_ratio(bbox, target_aspect_ratio, scale=1.0):
    """
    Expands a 2D bounding box to match a target aspect ratio (width / height),
    without shrinking. Optionally scales the final result uniformly, while preserving
    the target aspect ratio.

    Parameters:
    - bbox: np.ndarray, shape (2, 2)
        The input bounding box: [[xmin, ymin], [xmax, ymax]]
    - target_aspect_ratio: float
        Desired width / height ratio for the final bounding box.
    - scale: float, default=1.0
        A multiplicative factor to enlarge the final bounding box uniformly
        while maintaining the target aspect ratio. For example:
            scale = 1.0 → no change in size beyond aspect correction
            scale = 1.2 → 20% larger in both dimensions
            scale = 2.0 → double the size

    Returns:
    - expanded_bbox: np.ndarray, shape (2, 2)
        The adjusted bounding box with the desired aspect ratio and scale.
    """
    assert np.isscalar(target_aspect_ratio), "Aspect ratio must be scalar"
    assert target_aspect_ratio > 0, "Target aspect ratio must be positive"
    assert scale > 0, "Scale must be positive"

    bbox = np.array(bbox, dtype=float)
    assert bbox.shape == (2, 2), "bbox must be of shape (2, 2)"

    center = bbox.mean(axis=0)
    w, h = bbox[1] - bbox[0]
    current_aspect = w / h

    if current_aspect > target_aspect_ratio:
        # Too wide → expand height
        new_h = w / target_aspect_ratio
        new_w = w
    elif current_aspect < target_aspect_ratio:
        # Too tall → expand width
        new_w = h * target_aspect_ratio
        new_h = h
    else:
        new_w, new_h = w, h

    # Apply uniform scale while preserving aspect
    new_w *= scale
    new_h *= scale

    half_size = np.array([new_w / 2, new_h / 2])
    expanded_bbox = np.stack([center - half_size, center + half_size])

    return expanded_bbox



