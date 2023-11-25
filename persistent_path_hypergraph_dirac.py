


import numpy as np
import sympy




def common_member1(a, b):
    """
    Find common elements between two lists.

    Args:
        a (list): The first list.
        b (list): The second list.

    Returns:
        list: A list containing common elements between `a` and `b`.
    """
    result = [i for i in a if i in b]
    return result



## codes to remove elements of list from another, used below

def RemoveList(List, Subtract):
    """
    Remove elements from one list that are present in another list.

    Args:
        List (list): The list from which elements will be removed.
        Subtract (list): The list containing elements to be removed from `List`.

    Returns:
        None
    """
    for i in List[:]:
      if i in Subtract:
        List.remove(i)


# def Remove0(lst):
#     """
#     Remove duplicates from a list while preserving the order.

#     Args:
#         lst (list): The input list.

#     Returns:
#         list: A list with duplicate elements removed.
#     """
#     buffer = []
#     for i in lst:
#       buffer = buffer + [i]
#     return ([buffer(i) for i in {*[tuple(i) for i in buffer]}])


def Remove_Rep(lst):
    """
    Remove duplicate elements from a list.

    Args:
        lst (list): The input list.

    Returns:
        list: A list with duplicate elements removed.
    """
    return [*set(lst)]



def Flat_List(nestedlist):
    """
    Flatten a nested list.

    Args:
        nestedlist (list): The nested list to be flattened.

    Returns:
        list: A flat list containing all elements from the nested list.
    """
    flatlist=[]
    for sublist in nestedlist:
        for element in sublist:
            flatlist.append(element)
    return flatlist



def Remove_Repeated(k):
    """
    Remove repeated elements from a list while preserving the order.

    Args:
        k (list): The input list.

    Returns:
        list: A list with repeated elements removed.
    """
    new_k = []
    for elem in k:
        if elem not in new_k:
            new_k.append(elem)
    return new_k



# computing the boundary matrix
# ImageBasis are converted to IntegerLabel form, and the calculations are carried in the given order.
# The output is actually the TRANSPOSE OF THE BOUDNARY MATRIX


def I_thElementDropping(a):
    """
    Remove the i-th element from a list.

    Args:
        a (list): The input list.

    Returns:
        list: A list with the i-th element removed.
    """
    cartesian_length = len(a)
    A = []
    for i in range(cartesian_length):
      B = []
      B = a[:i] + a[i+1:] # a[:i] drop the element in the i-th position
      A.extend([B])
    return A




def Remove0(lst):
    """
    Remove duplicate lists from the input list 'lst' while preserving the original order.

    Args:
        lst (list): The input list containing lists.

    Returns:
        list: A new list with duplicate lists removed, preserving the original order.
    """
    return ([list(i) for i in {*[tuple(i) for i in lst]}])




def BoundaryMapOnFiniteSetsListLabeling(SourceBasis, ImageBasis):
    """
    Compute the boundary matrix.

    Args:
        SourceBasis (list of lists): List of source space basis elements.
        ImageBasis (list): List of image space basis elements.

    Returns:
        list of lists: The the boundary matrix.
    """
    BoundaryMatrix = [] # this is actually the transpose since we will use append method,
                        #remember the matrix is constructed by appending columns, but python append rows better
    SourceDim = len(SourceBasis)
    cartesian_length = len(SourceBasis[0]) # in what V^m SourceBasis is subset of its standard basis, assuming all in V^m
    ImageDim = len(ImageBasis)
    BufferVector = [0 for i in range(len(ImageBasis))]
    for i in range(SourceDim):
      BufferVector = [t * 0 for t in BufferVector] # trick reference https://stackoverflow.com/questions/35166633/how-do-i-multiply-each-element-in-a-list-by-a-number
      # compute image of an element in the source basis, and pick the element in the imagebasis with nonzero coefficients
      ImageOfBaseElementNonZeroComponents = I_thElementDropping(SourceBasis[i]) # it has the same lenght as SourceBasis[0], or generally SourceBasis elements
      # convert the above subcollection of basis elements to IntegerLabel, imaging V^m as digits in base |V|= size of V
      #### removed: ImageOfBaseElementConvertedToIntegerLables = ListBaseToNumber(ImageOfBaseElementNonZeroComponents, V_size)
      # find the corresponding index of the above elements in the ImageBasis list and assign the coefficient in BufferVector, so BufferVector is the coefficient matrix
      for j in range(cartesian_length):
        BufferVector[ImageBasis.index(ImageOfBaseElementNonZeroComponents[j])] = BufferVector[ImageBasis.index(ImageOfBaseElementNonZeroComponents[j])] + (-1)**j # Originally above, we used ImageOfBaseElementConvertedToIntegerLables
      BoundaryMatrix.extend([BufferVector])
    return np.array(BoundaryMatrix).T.tolist()




def FromSourceBasisGenerateImageSmallestSubbasisListLabeling(SourceBasis):
    """
    Generate the smallest subset of standard basis elements containing the image of a given source subbasis.

    Args:
        SourceBasis (list of lists): List of source basis elements.

    Returns:
        list: The smallest subset of standard basis elements containing the image of the source subbasis.
    """
    ImageBasis = []
    SourceDim = len(SourceBasis)
    for i in range(SourceDim):
      a = []
      # compute basis element image and pick basis elements with nonzero coefficients
      a = I_thElementDropping(SourceBasis[i])
      # convert to IntegerLabel form
      ### remove a = ListBaseToNumber(a, V_size)
      # append, though there might be some repeated values
      ImageBasis.extend(a)
    # remove repeated using numpy
    b = Remove0(ImageBasis)
    return b




def ProjectionBoundaryMapOnFiniteSetsListLabeling(SourceBasis, ImageBasis):
    """
    Compute a projection to the space generated by ImageBasis.

    Args:
        SourceBasis (list of lists): List of source basis elements.
        ImageBasis (list): List of image basis elements.

    Returns:
        list of lists: The the boundary matrix for the projection map.
    """
    BoundaryMatrix = [] # this is actually the transpose since we will use append method,
                        #remember the matrix is constructed by appending columns, but python append rows better
    SourceDim = len(SourceBasis)
    cartesian_length = len(SourceBasis[0]) # in what V^m SourceBasis is subset of its standard basis, assuming all in V^m
    ImageDim = len(ImageBasis)
    common_elements = []
    BufferVector = [0 for i in range(len(ImageBasis))]
    for i in range(SourceDim):
      BufferVector = [t * 0 for t in BufferVector] # trick reference https://stackoverflow.com/questions/35166633/how-do-i-multiply-each-element-in-a-list-by-a-number
      # compute image of an element in the source basis, and pick the element in the imagebasis with nonzero coefficients
      ImageOfBaseElementNonZeroComponents = I_thElementDropping(SourceBasis[i]) # it has the same lenght as SourceBasis[0], or generally SourceBasis elements
      # convert the above subcollection of basis elements to IntegerLabel, imaging V^m as digits in base |V|= size of V
      #### removed: ImageOfBaseElementConvertedToIntegerLables = ListBaseToNumber(ImageOfBaseElementNonZeroComponents, V_size)
      # find the corresponding index of the above elements in the ImageBasis list and assign the coefficient in BufferVector, so BufferVector is the coefficient matrix
      common_elements = common_member1(ImageOfBaseElementNonZeroComponents, ImageBasis)
      for base_element in common_elements:
        BufferVector[ImageBasis.index(base_element)] = BufferVector[ImageBasis.index(base_element)] + (-1)**ImageOfBaseElementNonZeroComponents.index(base_element) # Originally above, we used ImageOfBaseElementConvertedToIntegerLables
      BoundaryMatrix.extend([BufferVector])
    return np.array(BoundaryMatrix).T.tolist()



### find common elements of two python lists

def common_member1(a, b):
  result = [i for i in a if i in b]
  return result

### creater 2d zero python list

def PythonZeros(rows, columns):
    """
    Create a 2D list filled with zeros.

    Args:
        rows (int): Number of rows in the list.
        columns (int): Number of columns in the list.

    Returns:
        list of lists: A 2D list filled with zeros.
    """
    return [[0 for i in range(columns)] for j in range(rows)]


def Connected_to(graph, vertex): 
    """
    Find vertices connected to a given vertex in a graph.

    Args:
        graph (list of lists): The graph given as list of edges, edges are list of vertices.
        vertex: The vertex to find connections for.

    Returns:
        list: A list of vertices connected to the input vertex.
    """
    con_v = []
    for edge in graph:
      if vertex in edge:
        con_v = con_v + edge
    return (con_v)



## Same as above two codes but to return vertices only

def Ordered_DFS_VRegular_Paths_ListForm_Incidence_Vertex(graph, start, Length, paths = [], V=[], E=[]):
    """
    Perform a depth-first search for finding V-regular paths in a graph.
    V-regular = vertex-regular = anchor sequences where v_i not equal to v_i+1

    Args:
        graph (list of lists): The graph is given as list of edges, which are list of vertices.
        start: The starting vertex for the search.
        Length (int): The desired length of the paths.
        paths (list of lists): Accumulator for storing the found paths.
        V (list): List of vertices visited so far.
        E (list): List of edges visited so far.

    Returns:
        list of lists: List of V-regular paths found in the graph.
    """
    V1 = V + [start]
    if len(V1) > len(paths):
      paths = paths + [[]]
    paths[len(V1)-1] = paths[len(V1)-1] + [V1]
    if len(V1) == Length:
      return paths
    size = len(graph)
    for i in range(size):
      buffer = []
      buffer = buffer + graph[i]
      if start in buffer:
        E1 = E + [i]
        RemoveList(buffer,[start])
        for s in buffer:
          paths = Ordered_DFS_VRegular_Paths_ListForm_Incidence_Vertex(graph, s, Length, paths, V1, E1)
    return paths




#####################################
#####################################
#####################################
#####################################
#####################################




def VRegular_Paths_Vertex(Graph, VertexSet, Length):
    """
    Find V-regular paths starting from vertices in the given VertexSet.

    Args:
        Graph (dict): Dictionary representation of the graph with vertices as keys and adjacent vertices as values.
        VertexSet (list): List of vertices to start from.
        Length (int): Length of the V-regular paths to search for.

    Returns:
        list of lists: List of V-regular paths of the specified length.
    """
    graph_paths = [[] for i in range(Length)] # []
    for vertex in VertexSet:
      path_vertex = Ordered_DFS_VRegular_Paths_ListForm_Incidence_Vertex(Graph, vertex, Length)
      # if len(path_vertex) > len(graph_paths):
      #   for i in range(len(path_vertex) - len(graph_paths)):
      #     graph_paths = graph_paths + [[]]
      for i in range(len(path_vertex)):
        graph_paths[i] = graph_paths[i] + Remove_Repeated(path_vertex[i])
    return graph_paths

# Ordered lists of paths by length of lengths less than or equal to provided from start

def Ordered_DFS_VRegular_Paths_Vertex_Digraph(graph, start, Length, paths=[], V = []):
    """
      Perform depth-first search for finding V-regular paths in a directed graph.

      Args:
          graph (dict): Dictionary representation of the graph with vertices as keys and lists of outgoing vertices as values.
          start: The starting vertex for the search.
          Length (int): The desired length of the paths.
          paths (list of lists): Accumulator for storing the found paths.
          V (list): List of vertices visited so far.

      Returns:
          list of lists: List of V-regular paths found in the directed graph.
      """
    V1 = V + [start] # add vertex to the current sequence
    if len(V1) > len(paths):
      paths = paths + [[]]
    paths[len(V1)-1] = paths[len(V1)-1] + [V1]
    if len(V1) == Length:
      return paths
    Buffer = []
    Buffer = Buffer + graph[start]
    RemoveList(Buffer, [start]) # make regular
    for new_vertex in Buffer:
      paths = Ordered_DFS_VRegular_Paths_Vertex_Digraph(graph, new_vertex, Length, paths, V1)
    return paths



def VRegular_Paths_Vertex_Digraph(Graph, VertexSet, Length):
    """
      Find V-regular paths in a directed graph starting from vertices in the given VertexSet.

      Args:
          Graph (dict): Dictionary representation of the directed graph with vertices as keys and lists of outgoing vertices as values.
          VertexSet (list): List of vertices to start from.
          Length (int): Length of the V-regular paths to search for.

      Returns:
          list of lists: List of V-regular paths of the specified length in the directed graph.
      """
    graph_paths = [[] for i in range(Length)] # []
    for vertex in VertexSet:
      path_vertex = Ordered_DFS_VRegular_Paths_Vertex_Digraph(Graph, vertex, Length)
      # if len(path_vertex) > len(graph_paths):
      #   for i in range(len(path_vertex) - len(graph_paths)):
      #     graph_paths = graph_paths + [[]]
      for i in range(len(path_vertex)):
        graph_paths[i] = graph_paths[i] + Remove_Repeated(path_vertex[i])
    return graph_paths




### return span of intersection of two spaces given their spans in columns of matrices

def Intersection_Of_Spaces(ASpan, BSpan): 
    """
      Compute the basis of the intersection of two vector spaces given their generators.

      Args:
          ASpan (matrix of generators, column form): Basis of the first vector space.
          BSpan (same): Basis of the second vector space.

      Returns:
          matrix: Basis of the intersection of the two vector spaces in the first space.
      """
    # return basis of intersection of two spaces ASpan and BSpan, where these two are generating vectors
    # The return basis resides in ASpan space
    A_matrix = sympy.Matrix(ASpan)
    Acols = A_matrix.cols
    B_matrix = sympy.Matrix(BSpan)
    Intersection_Matrix = A_matrix.row_join(-B_matrix)
    NullAuxiliary = Intersection_Matrix.nullspace()
    ColumnsOutOfNullBasis = sympy.Matrix([[NullAuxiliary[i][j] for j in range(Acols)] for i in range(len(NullAuxiliary))]) # span of nullspace of A
    BasisAuxiliary = ColumnsOutOfNullBasis.rref()[0]
    AuxiliaryRank = BasisAuxiliary.rank()
    if AuxiliaryRank != 0:
      return BasisAuxiliary.T[:,0:AuxiliaryRank]
    return sympy.zeros(Acols, 1)


def Intersection_Of_Spaces2(ASpan, BSpan): # less memory used, but hard to read
    """
      Compute the basis of the intersection of two vector spaces given their basis (alternative method).

      Args:
          ASpan (matrix of generators, column form): Basis of the first vector space.
          BSpan (same): Basis of the second vector space.

      Returns:
          matrix: Basis of the intersection of the two vector spaces.
      """
    A_matrix = sympy.Matrix(ASpan)
    Acols = A_matrix.cols
    B_matrix = sympy.Matrix(BSpan)
    Intersection_Matrix = A_matrix.row_join(-B_matrix)
    BasisAuxiliary = sympy.Matrix([[Intersection_Matrix.nullspace()[i][j] for j in range(Acols)] for i in range(len(Intersection_Matrix.nullspace()))]).rref()[0]
    AuxiliaryRank = BasisAuxiliary.rank()
    if AuxiliaryRank != 0:
      return BasisAuxiliary.T[:,0:AuxiliaryRank]
    return sympy.zeros(Acols, 1)

###########################
###########################

# fast ones, numerical

def Fast_Intersection_Of_Spaces(ASpan, BSpan): # coefficient basis of the span of ASpan
    """
      Compute the basis of the intersection of two vector spaces given their basis (numerically optimized).

      Args:
          ASpan (matrix of generators, column form): Basis of the first vector space.
          BSpan (same): Basis of the second vector space.

      Returns:
          matrix: Basis of the intersection of the two vector spaces.
      """
    A_matrix = np.array(ASpan)
    Acols = A_matrix.shape[len(A_matrix.shape)-1]
    B_matrix = np.array(BSpan)
    Intersection_Matrix = np.concatenate((A_matrix, -B_matrix), axis=1)
    NullAuxiliary = sp.linalg.null_space(Intersection_Matrix, rcond = 0.00000000000000000000000000000000001)
    if NullAuxiliary.shape[1] != 0:
      ColumnsOutOfNullBasis = NullAuxiliary[0:Acols,:] # span of nullspace of A
      q, r = sp.linalg.qr(ColumnsOutOfNullBasis)
      if q.shape[1] != 0:
        q[np.abs(q) < 0.0000000000001] = 0 # tolerance to make it zero
        AuxiliaryRank = np.linalg.matrix_rank(q) # though q.shape returns rank in general, it does not if q is zero matrix
        if AuxiliaryRank != 0:
          return q # the basis
      return np.zeros([Acols, 1])
    return np.zeros([Acols, 1])



###########################
###########################
###########################
###########################


# return the matrix representaion from D=T^-1(B) -> B given T and the embdedding matrix of B, and coefficient basis of D

def Matrix_Representation_of_Inverse_Image_and_Basis(T, B):
    """
      Compute the matrix representation of the inverse image and the basis of the inverse image.

      Args:
          T (matrix): Matrix representation of the linear operator T.
          B (same): Target space spanning matrix.

      Returns:
          list: A list containing the basis of the inverse image and the matrix representation of the restriction of the operator on the inverse image with this basis.
      """
    # T is the matrix representation of the linear operator T
    # B is the target space spanning matrix
    # return the basis of the inverse image M in the source, 
    # Aux.inv is the representation of the restriction of the operator on M with this basis
    T_matrix = sympy.Matrix(T)
    B_matrix = sympy.Matrix(B)
    M = Intersection_Of_Spaces2(T_matrix, B_matrix) # D Coefficient Basis
    Aux = B_matrix.T*B_matrix
    return [M, Aux.inv(method="LU")*B_matrix.T*T_matrix*M]

def Matrix_Representation_of_Inverse_Image(T, B):
    """
      Compute the matrix representation of the inverse image.

      Args:
          T (list of lists): Matrix representation of the linear operator T.
          B (list of lists): Target space spanning matrix.

      Returns:
          list: The matrix representation of the restriction of the operator on the inverse image with a basis.
      """
    # T is the matrix representation of the linear operator T
    # B is the target space spanning matrix
    # return the basis of the inverse image M in the source, 
    # Aux.inv is the representation of the restriction of the operator on M with this basis
    T_matrix = sympy.Matrix(T)
    B_matrix = sympy.Matrix(B)
    M = Intersection_Of_Spaces2(T_matrix, B_matrix) # D Coefficient Basis
    Aux = B_matrix.T*B_matrix
    return Aux.inv(method="LU")*B_matrix.T*T_matrix*M

def Basis_Representation_of_Inverse_Image(T, B):
    """
      Compute the basis of the inverse image.

      Args:
          T (list of lists): Matrix representation of the linear operator T.
          B (list of lists): Target space spanning matrix.

      Returns:
          list: The basis of the inverse image in the source space.
      """
    # T is the matrix representation of the linear operator T
    # B is the target space spanning matrix
    # return the basis of the inverse image M in the source, 
    # Aux.inv is the representation of the restriction of the operator on M with this basis
    T_matrix = sympy.Matrix(T)
    B_matrix = sympy.Matrix(B)
    M = Intersection_Of_Spaces2(T_matrix, B_matrix) # D Coefficient Basis
    return M


###########################
###########################

# fast ones, numerical

def Fast_Matrix_Representation_of_Inverse_Image_and_Basis(T, B):
    """
      Compute the matrix representation of the inverse image and the basis of the inverse image (numerically optimized).

      Args:
          T (list of lists): Matrix representation of the linear operator T.
          B (list of lists): Target space spanning matrix.

      Returns:
          list: A list containing the basis of the inverse image and the matrix representation of the restriction of the operator on the inverse image with this basis.
      """
    T_matrix = np.array(T)
    B_matrix = np.array(B)
    M = Fast_Intersection_Of_Spaces(T_matrix, B_matrix) # D Coefficient Basis
    Aux = B_matrix.T@B_matrix
    return [M, np.linalg.inv(Aux)@B_matrix.T@T_matrix@M]

def Fast_Matrix_Representation_of_Inverse_Image(T, B):
    """
      Compute the matrix representation of the inverse image (numerically optimized).

      Args:
          T (list of lists): Matrix representation of the linear operator T.
          B (list of lists): Target space spanning matrix.

      Returns:
          list: The matrix representation of the restriction of the operator on the inverse image with a basis.
      """
    T_matrix = np.array(T)
    B_matrix = np.array(B)
    M = Fast_Intersection_Of_Spaces(T_matrix, B_matrix) # D Coefficient Basis
    Aux = B_matrix.T@B_matrix
    return np.linalg.inv(Aux)@B_matrix.T@T_matrix@M

def Fast_Basis_Representation_of_Inverse_Image(T, B):
    """
      Compute the basis of the inverse image (numerically optimized).

      Args:
          T (list of lists): Matrix representation of the linear operator T.
          B (list of lists): Target space spanning matrix.

      Returns:
          list: The basis of the inverse image in the source space.
      """
    T_matrix = np.array(T)
    B_matrix = np.array(B)
    M = Fast_Intersection_Of_Spaces(T_matrix, B_matrix) # D Coefficient Basis
    return M


###########################
###########################
###########################
###########################

## embedding matrix of subspace with basis extended to Basis
def Embedding_Matrix(Subbasis, Basis):
    """
      Compute the embedding matrix of a subspace with its basis extended to a larger basis.

      Args:
          Subbasis (list): Basis of the subspace.
          Basis (list): Larger basis that includes the subspace basis.

      Returns:
          list of lists: Embedding matrix that represents the embedding of the subspace into the larger space.
      """
    # given two labels of basis, one is sub-collection of the other
    # return a matrix representation of the embedding, namely 0's 1's kind of matrix
    SizeSubbasis = len(Subbasis)
    SizeBasis = len(Basis)
    EmbeddingMatrix = np.zeros((SizeBasis, SizeSubbasis)).tolist() # buffer of zeros, prefer to use np to avoid row cases [[]] == [] issue when dealing with others
    for i in Subbasis:
        EmbeddingMatrix[Basis.index(i)][Subbasis.index(i)] = 1
    return EmbeddingMatrix


## embedding matrix of subspace with basis extended to Basis
def Embedding_Matrix2(Subbasis, Basis):
    """
      Compute the embedding matrix of a subspace with its basis extended to a larger basis (alternative method).

      Args:
          Subbasis (list): Basis of the subspace.
          Basis (list): Larger basis that includes the subspace basis.

      Returns:
          list of lists: Embedding matrix that represents the embedding of the subspace into the larger space.
      """
    # SizeSubbasis = len(Subbasis)
    SizeBasis = len(Basis)
    Common = []
    Common = common_member1(Subbasis, Basis)
    Common_Size = len(Common)
    EmbeddingMatrix = np.zeros((SizeBasis, Common_Size)).tolist() # buffer of zeros, prefer to use np to avoid row cases [[]] == [] issue when dealing with others
    for i in Common:
      # if i in Basis:
      EmbeddingMatrix[Basis.index(i)][Common.index(i)] = 1
    return EmbeddingMatrix




def CommonSpaceBasis(Basis1, Basis2):
    """
      Find the intersection of two sets of basis.

      Args:
          Basis1 (list): First set of basis.
          Basis2 (list): Second set of basis.

      Returns:
          list: List of basis elements that are common to both sets.
      """
    # given two labels of basis, find the intersection of these two
    Common = []
    Common.extend(Basis1)
    Common.extend(Basis2)
    Common = Remove0(Common)
    return Common





# Given sequence of submodules (A_n) of a complex (Omega_n) such that A_n submodule of Omega_n
# return the largest subcomplex contained in the this sequence.

def SubbSpaces_Complex_From_Generating_Basis(SpanBasis):
    """
      Compute the largest subcomplex contained in a sequence of submodules.

      Args:
          SpanBasis (list of lists): List of submodules' generating bases.

      Returns:
          list: A list containing the anchor sequences basis, path complex component, and boundary maps.
      """
    An = SpanBasis
    Length = len(An)
    Mn = [[] for i in range(Length)] # Dn+1 -> An+1
    PCn = [[] for i in range(Length)]

    # initial case
    if An[0] == []:
      Mn[0] = sympy.Matrix([0])
      PCn[0] = sympy.Matrix([0])
    else:
      Mn[0] = sympy.eye(len(An[0])) # D0 = A0 always, a choice made
      PCn[0] = sympy.Matrix([[0 for i in range(len(An[0]))]])

    # rest
    for i in range(Length - 1):
      EmAnCommon = []
      BImAn1 = []
      CommonSpace = []
      ImAn1 = []
      Buffer = []
      if An[i+1] == []:
        PCn[i+1] = sympy.zeros(Mn[i].cols, 1)
        if An[i] == []:
          Mn[i+1] = sympy.Matrix([0])
          EmAnCommon = sympy.Matrix([0])
          BImAn1 = sympy.Matrix([0])
        else:
          Mn[i+1] = sympy.Matrix([0])
          EmAnCommon = sympy.eye(len(An[i]))
          BImAn1 = sympy.zeros(len(An[i]), 1)
      else:
        ImAn1 = FromSourceBasisGenerateImageSmallestSubbasisListLabeling(An[i+1])
        CommonSpace = CommonSpaceBasis(An[i], ImAn1)
        BImAn1 = sympy.Matrix(BoundaryMapOnFiniteSetsListLabeling(An[i+1], CommonSpace))
        if An[i] == []:
          EmAnCommon = sympy.zeros(len(CommonSpace), 1)
          Mn[i+1] = sympy.Matrix(Basis_Representation_of_Inverse_Image(BImAn1 , EmAnCommon*Mn[i]))
          if sympy.det(Mn[i+1].T*Mn[i+1]) != 0:
            Mn[i+1] = Mn[i+1].QRdecomposition()[0]
          PCn[i+1] = sympy.zeros(Mn[i].cols, Mn[i+1].cols)
        else:
          EmAnCommon = sympy.Matrix(Embedding_Matrix(An[i], CommonSpace))
          if sympy.det(Mn[i].T*Mn[i]) != 0:
            Buffer = Matrix_Representation_of_Inverse_Image_and_Basis(BImAn1 , EmAnCommon*Mn[i])
            if sympy.det(Buffer[0].T*Buffer[0]) != 0:
              Q, R = sympy.Matrix(Buffer[0]).QRdecomposition()
              Mn[i+1] = Q
              PCn[i+1] = sympy.Matrix(Buffer[1])*R.inv(method="LU")
            else:
              Mn[i+1] = sympy.Matrix(Buffer[0])
              PCn[i+1] = sympy.Matrix(Buffer[1])
          else:
            Mn[i+1] = sympy.Matrix(Basis_Representation_of_Inverse_Image(BImAn1 , EmAnCommon*Mn[i]))
            if sympy.det(Mn[i+1].T*Mn[i+1]) != 0:
              Mn[i+1] = Mn[i+1].QRdecomposition()[0]
            PCn[i+1] = sympy.zeros(Mn[i].cols, Mn[i+1].cols) # Nothing will be weird since Ans are not empty things, Mn will be in worse case sequence of zeros reflecting An dimension



      # ####### Test #######
      # ####### Test #######
      # ####### Test #######
      # ####### Test #######

      # print("Commutativity test")
      # print(BImAn1*Mn[i+1] == EmAnCommon*Mn[i]*PCn[i+1])
      # TestPCn = PCn[i]*PCn[i+1]
      # print("Complex test")
      # print(TestPCn + TestPCn == TestPCn)

      # ####### Test #######
      # ####### Test #######
      # ####### Test #######
      # ####### Test #######


    return [An, Mn, PCn] # An anchor sequences basis, Mn path complex component, PCn boundary maps


###################################################
###################################################
###################################################
###################################################

def Persistence_SubbSpaces_Complex_From_Generating_Basis(SpanBasis1, SpanBasis2): # like not all the needed maps are really constructed between the complexes
    """
    Compute the persistence subcomplexes from generating bases of two complexes.

    Args:
        SpanBasis1 (list of lists): Generating bases of the first complex.
        SpanBasis2 (list of lists): Generating bases of the second complex.

    Returns:
        list: A list containing various components of the persistence subcomplexes.
        
        # An, Bn anchor sequences
        # PCAn, PCBn are the boundary maps of A and B complexes,
        # MAn, MBn path complexes,
        # Vn embedding maps of auxiliary complex to compute persistent dirac,
        # SCn boundary maps of the auxiliary complex
        # Tn the map C_n+1 -> A_n, where C_n+1 is the n+1 component of the auxiliary complex
    """
    [An, MAn, PCAn] = SubbSpaces_Complex_From_Generating_Basis(SpanBasis1)
    [Bn, MBn, PCBn] = SubbSpaces_Complex_From_Generating_Basis(SpanBasis2)
    Length = len(SpanBasis1)
    Vn = [[] for i in range(Length)] # Cn -> DBn
    Tn = [[] for i in range(Length)] # the diagonal maps Cn - > DAn

    # initial case
    Vn[0] = MBn[0]
    Tn[0] = PCBn[0]


    for i in range(Length - 1):
      EmAnCommon = []
      BImBn1 = []
      CommonSpace = []
      ImBn1 = []
      Buffer = []
      if Bn[i+1] == []:
        Vn[i+1] = sympy.Matrix([0])
        Tn[i+1] = sympy.zeros(MAn[i].cols, Vn[i+1].cols)
        if An[i] == []:
          EmAnCommon = sympy.Matrix([0])
          BImBn1 = sympy.Matrix([0])
        else:
          EmAnCommon = sympy.eye(len(An[i]))
          BImBn1 = sympy.zeros(len(An[i]), 1)
      else:
        ImBn1 = FromSourceBasisGenerateImageSmallestSubbasisListLabeling(Bn[i+1])
        CommonSpace = CommonSpaceBasis(An[i], ImBn1)
        BImBn1 = sympy.Matrix(BoundaryMapOnFiniteSetsListLabeling(Bn[i+1], CommonSpace))
        if An[i] == []:
          EmAnCommon = sympy.zeros(len(CommonSpace), 1)
          if sympy.det(MBn[i+1].T*MBn[i+1]) != 0:
              Vn[i+1] = Basis_Representation_of_Inverse_Image(BImBn1*MBn[i+1], sympy.zeros(BImBn1.rows,1))
              if sympy.det(Vn[i+1].T*Vn[i+1]) != 0:
                Vn[i+1] = Vn[i+1].QRdecomposition()[0]
              Tn[i+1] = sympy.zeros(1, Vn[i+1].cols)
          else: # MBn = 0 space
            Vn[i+1] = sympy.Matrix([0])
            Tn[i+1] = sympy.zeros(MAn[i].cols, Vn[i+1].cols)
        else:
          EmAnCommon = sympy.Matrix(Embedding_Matrix(An[i], CommonSpace))
          X_Buffer = []
          X_Buffer = EmAnCommon*MAn[i]
          if sympy.det(MAn[i].T*MAn[i]) != 0:
            Buffer = Matrix_Representation_of_Inverse_Image_and_Basis(BImBn1*MBn[i+1], X_Buffer)
            # Vn[i+1] = sympy.Matrix(Buffer[0])
            # Tn[i+1] = sympy.Matrix(Buffer[1])
            # Buffer = Matrix_Representation_of_Inverse_Image_and_Basis(BImAn1 , EmAnCommon*Mn[i])
            if sympy.det(Buffer[0].T*Buffer[0]) != 0:
              Q, R = sympy.Matrix(Buffer[0]).QRdecomposition()
              Vn[i+1] = Q
              Tn[i+1] = sympy.Matrix(Buffer[1])*R.inv(method="LU")
            else:
              Vn[i+1] = sympy.Matrix(Buffer[0])
              Tn[i+1] = sympy.Matrix(Buffer[1])
          else: # MAn = 0 space
            if sympy.det(MBn[i+1].T*MBn[i+1]) != 0:
              Vn[i+1] = sympy.Matrix(Basis_Representation_of_Inverse_Image(BImBn1*MBn[i+1], X_Buffer))
              if sympy.det(Vn[i+1].T*Vn[i+1]) != 0:
                Vn[i+1] = Vn[i+1].QRdecomposition()[0]
              Tn[i+1] = sympy.zeros(MAn[i].cols, Vn[i+1].cols)
            else: # MBn = 0 space
              Vn[i+1] = sympy.Matrix([0])
              Tn[i+1] = sympy.zeros(MAn[i].cols, Vn[i+1].cols)



      # ####### Test #######
      # ####### Test #######
      # ####### Test #######
      # ####### Test #######

      # print("Commutativity test")
      # print(BImBn1*MBn[i+1]*Vn[i+1] == EmAnCommon*MAn[i]*Tn[i+1])

      # ####### Test #######
      # ####### Test #######
      # ####### Test #######
      # ####### Test #######

    # initial
    SCn = [[] for i in range(Length)] # the boundary maps
    SCn[0] = PCBn[0]

    # rest
    for i in range(Length - 1):
      if sympy.det(Vn[i].T*Vn[i]) != 0:
        SCn[i+1] = sympy.Matrix(Matrix_Representation_of_Inverse_Image(PCBn[i+1]*Vn[i+1], Vn[i]))
      else:
        SCn[i+1] = sympy.zeros(Vn[i].cols, Vn[i+1].cols) # hopefully this will fix the zero spaces

      # ####### Test #######
      # ####### Test #######
      # ####### Test #######
      # ####### Test #######

      # print("Commutativity test")
      # print(Vn[i]*SCn[i+1] == PCBn[i+1]*Vn[i+1])
      # TestPCn = SCn[i]*SCn[i+1]
      # print("Complex test")
      # print(TestPCn + TestPCn == TestPCn)

      # ####### Test #######
      # ####### Test #######
      # ####### Test #######
      # ####### Test #######


    return [PCAn, PCBn, SCn, Tn, An, MAn, Bn, MBn, Vn] 
    # An, Bn anchor sequences
    # PCAn, PCBn are the boundary maps of A and B complexes,
    # MAn, MBn path complexes,
    # Vn embedding maps of auxiliary complex to compute persistent dirac,
    # SCn boundary maps of the auxiliary complex
    # Tn the map C_n+1 -> A_n, where C_n+1 is the n+1 component of the auxiliary complex 

###################################################
###################################################
###################################################
###################################################




# Hypergraphs

def Path_Complex_Hypergraph(graph, vertex_set, Length):
    """
      Compute the path complex of a hypergraph.

      Args:
          graph (list of lists): Hypergraph represented as a list of hyperedges.
          vertex_set (list): List of vertices in the hypergraph.
          Length (int): Length of the desired path complex.

      Returns:
          list: A list containing various components of the path complex.
      """
    ### first compute the basis (regular paths) An
    SpanBasis = VRegular_Paths_Vertex(graph, vertex_set, Length)
    return SubbSpaces_Complex_From_Generating_Basis(SpanBasis)


def Persistence_Path_Complex_Hypergraph(graph1, vertex_set1, graph2, vertex_set2, Length):
    """
      Compute the persistence path complex of two hypergraphs.

      Args:
          graph1 (list of lists): First hypergraph represented as a list of hyperedges.
          vertex_set1 (list): List of vertices in the first hypergraph.
          graph2 (list of lists): Second hypergraph represented as a list of hyperedges.
          vertex_set2 (list): List of vertices in the second hypergraph.
          Length (int): Length of the desired path complex.

      Returns:
          list: A list containing various components of the persistence path complex.
      """

    ### first compute the basis (regular paths) An
    SpanBasis1 = VRegular_Paths_Vertex(graph1, vertex_set1, Length)
    SpanBasis2 = VRegular_Paths_Vertex(graph2, vertex_set2, Length)

    return Persistence_SubbSpaces_Complex_From_Generating_Basis(SpanBasis1, SpanBasis2)


###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################


# digraphs

def Path_Complex_Digraph(graph, VertexSet, Length):
    """
      Compute the path complex of a directed graph (digraph).

      Args:
          graph (dict): Directed graph represented as an adjacency list.
          VertexSet (list): List of vertices in the directed graph.
          Length (int): Length of the desired path complex.

      Returns:
          list: A list containing various components of the path complex.
      """
    ### first compute the basis (regular paths) An
    SpanBasis = VRegular_Paths_Vertex_Digraph(graph, VertexSet, Length)
    return SubbSpaces_Complex_From_Generating_Basis(SpanBasis)


def Persistence_Path_Complex_Digraph(graph1, VertexSet1, graph2, VertexSet2, Length): # Make sure the cosistency of labeling of both graphs graph1 -> graph2
    """
      Compute the persistence path complex of two directed graphs (digraphs).

      Args:
          graph1 (dict): First directed graph represented as an adjacency list.
          VertexSet1 (list): List of vertices in the first directed graph.
          graph2 (dict): Second directed graph represented as an adjacency list.
          VertexSet2 (list): List of vertices in the second directed graph.
          Length (int): Length of the desired path complex.

      Returns:
          list: A list containing various components of the persistence path complex.
      """
    ### first compute the basis (regular paths) An
    SpanBasis1 = VRegular_Paths_Vertex_Digraph(graph1, VertexSet1, Length)
    SpanBasis2 = VRegular_Paths_Vertex_Digraph(graph2, VertexSet2, Length)

    return Persistence_SubbSpaces_Complex_From_Generating_Basis(SpanBasis1, SpanBasis2)


###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################





def Matrix_Nullity(A):
    """
      Compute the nullity of a matrix.

      Args:
          A (list of lists): Input matrix represented as a list of lists.

      Returns:
          int: Nullity of the matrix.
      """
    return sympy.Matrix(A).cols - np.linalg.matrix_rank(np.array(A).astype(np.float64))

###################################################
###################################################

def EigenVals(A):
    """
      Compute the eigenvalues of a matrix.

      Args:
          A (list of lists): Input matrix represented as a list of lists.

      Returns:
          list: List of eigenvalues (rounded to 3 decimal places).
      """
    return np.linalg.eigvals(np.array(A).astype(np.float64)).round(3).real.tolist()

###################################################
###################################################

# def Fiedler_Vals_tol(A, tol): # A is a list of nonnegative eigenvalues
#   A_new = np.array(A).astype(np.float64) # make sure it is in np
#   A_new[np.abs(A_new) < tol ] = 0
#   return np.amin( A_new[A_new > 0])

def Fiedler_Vals(A): # A is a list of nonnegative eigenvalues
    """
      Compute the Fiedler value of a matrix.

      Args:
          A (list of numbers): List of non-negative eigenvalues.

      Returns:
          float: Fiedler value of the matrix.
      """
    A_new = np.array(A).astype(np.float64).real # make sure it is in np
    buffer = A_new[A_new > 0]
    if len(buffer)>0:
        return np.amin(buffer)
    return 0


###################################################
###################################################


def Maxs(A): # A is a list of nonnegative eigenvalues
    """
      Compute the maximum eigenvalue of a matrix.

      Args:
          A (list of numbers): List of non-negative eigenvalues.

      Returns:
          float: Maximum eigenvalue of the matrix.
      """
    A_new = np.array(A).astype(np.float64).real # make sure it is in np
    # return np.amax( A_new[A_new > 0])
    buffer = A_new[A_new > 0]
    if len(buffer)>0:
        return np.amax(buffer)
    return 0



###################################################
###################################################

def Positive_Mean(A): # A is a list of nonnegative eigenvalues
    """
      Compute the positive mean of a list of nonnegative eigenvalues.

      Args:
          A (list): List of nonnegative eigenvalues.

      Returns:
          float: Positive mean of the eigenvalues.
      """
    A_new = np.array(A).astype(np.float64).real # make sure it is in np
    # return np.mean( A_new[A_new > 0])
    buffer = A_new[A_new > 0]
    if len(buffer)>0:
        return np.mean(buffer)
    return 0


###################################################
###################################################


def Positive_Sum(A): # A is a list of nonnegative eigenvalues
    """
      Compute the positive sum of a list of nonnegative eigenvalues.

      Args:
          A (list): List of nonnegative eigenvalues.

      Returns:
          float: Positive sum of the eigenvalues.
      """
    A_new = np.array(A).astype(np.float64).real # make sure it is in np
    # return np.sum( A_new[A_new > 0])
    buffer = A_new[A_new > 0]
    if len(buffer)>0:
        return np.sum(buffer)
    return 0

###################################################
###################################################


def Positive_Generalized_Mean(A): # A is a list of nonnegative eigenvalues
    """
      Compute the positive generalized mean of a list of nonnegative eigenvalues.

      Args:
          A (list): List of nonnegative eigenvalues.

      Returns:
          float: Positive generalized mean of the eigenvalues.
      """
    mean_value = Positive_Mean(A)
    A_new = np.array(A).astype(np.float64).real # make sure it is in np
    A_new = A_new - mean_value
    return np.mean(np.abs(A_new))


###################################################
###################################################



def Dirac_Nullity(LapsNulls, DownLapsNulls):
    """
      Compute the Dirac nullity from Laplacian nullities and DownLaplacian nullities.

      Args:
          LapsNulls (list): List of Laplacian nullities.
          DownLapsNulls (list): List of DownLaplacian nullities.

      Returns:
          list: List of Dirac nullities.
      """
    Dirac = [[] for i in range(len(LapsNulls))]
    for i in range(len(LapsNulls) - 1):
      Dirac[i] = sum(LapsNulls[:i+1]) + DownLapsNulls[i+1]
    Dirac[len(LapsNulls) - 1] = sum(LapsNulls)
    return Dirac

def Dirac_NonNegEigenVals(LapsEigens, DownsEigens):
    """
      Compute the Dirac nonnegative eigenvalues from Laplacian eigenvalues and DownLaplacian eigenvalues.

      Args:
          LapsEigens (list): List of Laplacian eigenvalues.
          DownsEigens (list): List of DownLaplacian eigenvalues.

      Returns:
          list: List of Dirac nonnegative eigenvalues.
      """
    Dirac_NonNegEigens = [[] for i in range(len(LapsEigens))]
    for i in range(len(LapsEigens) - 1):
      Dirac_NonNegEigens[i] = np.sqrt(Flat_List(LapsEigens[:i+1]) + DownsEigens[i+1]).tolist()
    Dirac_NonNegEigens[len(LapsEigens) - 1] = np.sqrt(Flat_List(LapsEigens)).tolist()
    return Dirac_NonNegEigens


def Dirac_Fiedler_Vals(Dirac_NonNegEigens):
    """
      Compute the Dirac Fiedler values from Dirac nonnegative eigenvalues.

      Args:
          Dirac_NonNegEigens (list): List of Dirac nonnegative eigenvalues.

      Returns:
          list: List of Dirac Fiedler values.
      """
    return [Fiedler_Vals(Dirac_NonNegEigens[i]) for i in range(len(Dirac_NonNegEigens))]

def Dirac_Maxs(Dirac_NonNegEigens):
    """
      Compute the Dirac maximum eigenvalues from Dirac nonnegative eigenvalues.

      Args:
          Dirac_NonNegEigens (list): List of Dirac nonnegative eigenvalues.

      Returns:
          list: List of Dirac maximum eigenvalues.
      """
    return [Maxs(Dirac_NonNegEigens[i]) for i in range(len(Dirac_NonNegEigens))]

def Dirac_Positive_Mean(Dirac_NonNegEigens):
    """
      Compute the Dirac positive mean from Dirac nonnegative eigenvalues.

      Args:
          Dirac_NonNegEigens (list): List of Dirac nonnegative eigenvalues.

      Returns:
          list: List of Dirac positive mean values.
      """
    return [Positive_Mean(Dirac_NonNegEigens[i]) for i in range(len(Dirac_NonNegEigens))]

def Dirac_Positive_Sum(Dirac_NonNegEigens):
    """
      Compute the Dirac positive sum from Dirac nonnegative eigenvalues.

      Args:
          Dirac_NonNegEigens (list): List of Dirac nonnegative eigenvalues.

      Returns:
          list: List of Dirac positive sum values.
      """
    return [Positive_Sum(Dirac_NonNegEigens[i]) for i in range(len(Dirac_NonNegEigens))]

def Dirac_Positive_Generalized_Mean(Dirac_NonNegEigens):
    """
      Compute the Dirac positive generalized mean from Dirac nonnegative eigenvalues.

      Args:
          Dirac_NonNegEigens (list): List of Dirac nonnegative eigenvalues.

      Returns:
          list: List of Dirac positive generalized mean values.
      """
    return [Positive_Generalized_Mean(Dirac_NonNegEigens[i]) for i in range(len(Dirac_NonNegEigens))]




def Initial_Dirac(Boundary_Matrix, Embedding_of_Source, Embedding_of_Target):
    """
      Compute the initial Dirac from boundary matrix and embeddings.

      Args:
          Boundary_Matrix (list of lists): Boundary matrix.
          Embedding_of_Source (list of lists): Embedding of source complex.
          Embedding_of_Target (list of lists): Embedding of target complex.

      Returns:
          list of lists: Initial Dirac.
      """
    Bound = np.array(Boundary_Matrix).astype('float64')
    Is = np.array(Embedding_of_Source).astype('float64')
    It = np.array(Embedding_of_Target).astype('float64')
    A = np.concatenate((np.zeros([Bound.shape[0], Bound.shape[0]]), Bound), axis = 1) # [0 B]
    if np.abs(np.linalg.det(Is.T@Is)) < 0.001:
      B = np.zeros_like(Bound).T
    else:
      B = np.linalg.inv((Is.conjugate().T@Is))@Bound.conjugate().T@(It.conjugate().T@It) # adjoint operator
    C = np.concatenate((B, np.zeros([Bound.shape[1], Bound.shape[1]])), axis = 1) # [B* 0]
    D = np.concatenate((A, C)) # Dirac0
    return D


def Next_Dirac(Previous_Dirac, Next_Boundary, Embedding_of_Source, Embedding_of_Target):
    """
      Compute the next Dirac from previous Dirac, next boundary, and embeddings.

      Args:
          Previous_Dirac (list of lists): Previous Dirac.
          Next_Boundary (list of lists): Next boundary.
          Embedding_of_Source (list of lists): Embedding of source complex.
          Embedding_of_Target (list of lists): Embedding of target complex.

      Returns:
          list of lists: Next Dirac.
      """
    A = np.array(Previous_Dirac).astype('float64')
    B = np.array(Next_Boundary).astype('float64')
    Is = np.array(Embedding_of_Source).astype('float64')
    It = np.array(Embedding_of_Target).astype('float64')
    # 1
    C = np.concatenate(\
    (np.zeros([A.shape[0] - B.shape[0], B.shape[1]]),\
      B), axis=0) # [[0], [B]]
    # 2
    D = np.concatenate((A, C), axis = 1) # [A C]
    # 3
    # E = B.conjugate().T # adjoint operator
    if np.abs(np.linalg.det(Is.T@Is)) < 0.001:
      E = np.zeros_like(B).T
    else:
      E = np.linalg.inv((Is.conjugate().T@Is))@B.conjugate().T@(It.conjugate().T@It) # adjoint operator
    F = np.concatenate((np.zeros([E.shape[0], A.shape[1]-E.shape[1]]), E), axis = 1) # [0 E]
    # 4
    G = np.concatenate((F, np.zeros([E.shape[0], B.shape[1]])), axis = 1) # [F 0] = [0 E 0]
    # 5
    H = np.concatenate((D, G), axis = 0) # [[D], [G]]
    return H


def Complex_Diracs(Boundary_Maps, Matrix_Representation_Embbedings): # include the first zero map
    """
      Compute a list of Diracs from boundary maps and matrix representation embeddings.

      Args:
          Boundary_Maps (list of lists): List of boundary maps.
          Matrix_Representation_Embbedings (list of lists): List of matrix representation embeddings.

      Returns:
          list of lists: List of Diracs.
      """
    Length = len(Boundary_Maps)
    D = [[] for i in range(Length - 1)]
    D[0] = Initial_Dirac(Boundary_Maps[1], Matrix_Representation_Embbedings[1], Matrix_Representation_Embbedings[0])
    for i in range(1, Length - 1):
      D[i] = Next_Dirac(D[i-1], Boundary_Maps[i+1], Matrix_Representation_Embbedings[i+1], Matrix_Representation_Embbedings[i])
    return D



###################################################
###################################################

def MatrixFeatures(A):
    """
      Compute various features of a matrix.

      Args:
          A (list of lists): Input matrix represented as a list of lists.

      Returns:
          dict: A dictionary containing matrix features.
      """
    Feature = {}
    Feature.update({'Null': Matrix_Nullity(A)})
    Feature.update({'Eigens': EigenVals(A)})
    Feature.update({'Positive_Sum':Positive_Sum(Feature['Eigens'])})
    Feature.update({'Positive_Generalized_Mean':Positive_Generalized_Mean(Feature['Eigens'])})
    Feature.update({'Maxs':Maxs(Feature['Eigens'])})
    Feature.update({'Fiedler_Vals':Fiedler_Vals(Feature['Eigens'])})
    Feature.update({'Positive_Mean':Positive_Mean(Feature['Eigens'])})
    return Feature

def Dirac_Features(DownLapsNulls, LapsNulls, DownsEigens, LapsEigens):
  Feature = {}
  Feature.update({'Null' : Dirac_Nullity(LapsNulls, DownLapsNulls)})
  Feature.update({'Eigens' : Dirac_NonNegEigenVals(LapsEigens, DownsEigens)})
  Dirac_NonNegEigens = Feature['Eigens']
  Feature.update({'Positive_Sum' : Dirac_Positive_Sum(Dirac_NonNegEigens)})
  Feature.update({'Positive_Generalized_Mean' : Dirac_Positive_Generalized_Mean(Dirac_NonNegEigens)})
  Feature.update({'Maxs' : Dirac_Maxs(Dirac_NonNegEigens)})
  Feature.update({'Fiedler_Vals' : Dirac_Fiedler_Vals(Dirac_NonNegEigens)})
  Feature.update({'Positive_Mean' : Dirac_Positive_Mean(Dirac_NonNegEigens)})
  return Feature




# Related Matrices

def UpLaplacians_Of_Complex(Boundary_Maps, Matrix_Representation_Embbedings): # f: V -> W linear, then upper laplacian is at W, given by f*f^t
    """
    Compute the upper Laplacians of a complex.

    Args:
        Boundary_Maps (list of arrays): List of boundary maps.
        Matrix_Representation_Embbedings (list of arrays): Mn in the above codes.

    Returns:
        list of arrays: List of upper Laplacians corresponding to each boundary map.
    """
    A = Boundary_Maps
    Reps = Matrix_Representation_Embbedings
    L = len(A)
    UpLaps = [ 0 for i in range(L)]
    for i in range(L-1):
      B = np.array(A[i+1]).astype(np.float64)
      B_hat = Adjoint_Operator(B, Reps[i+1], Reps[i])
      UpLaps[i] = B@B_hat
    UpLaps[L-1] = np.array(sympy.zeros(sympy.Matrix(A[L-1]).cols)).astype('float64')
    return UpLaps
###################################################
###################################################


def DownLaplacians_Of_Complex(Boundary_Maps, Matrix_Representation_Embbedings): # f: V -> W linear, then upper laplacian is at W, given by f*f^t
    """
    Compute the lower Laplacians of a complex.

    Args:
        Boundary_Maps (list of arrays): List of boundary maps.
        Matrix_Representation_Embbedings (list of arrays): List of matrix representation embeddings.

    Returns:
        list of arrays: List of lower Laplacians corresponding to each boundary map.
    """
    A = Boundary_Maps
    Reps = Matrix_Representation_Embbedings
    L = len(A)
    DownLaps = [ 0 for i in range(L)]
    for i in range(1, L):
      B = np.array(A[i]).astype(np.float64)
      B_hat = Adjoint_Operator(B, Reps[i], Reps[i-1])
      DownLaps[i] = B_hat@B
    DownLaps[0] = np.array(sympy.zeros(sympy.Matrix(A[1]).rows)).astype('float64')
    return DownLaps

###################################################
###################################################


def Laplacians_Of_Complex(Boundary_Maps, Matrix_Representation_Embbedings):
    """
    Compute the Laplacians of a complex as the sum of upper and lower Laplacians.

    Args:
        Boundary_Maps (list of arrays): List of boundary maps.
        Matrix_Representation_Embbedings (list of arrays): List of matrix representation embeddings.

    Returns:
        list of arrays: List of Laplacians corresponding to each boundary map.
    """
    Up = UpLaplacians_Of_Complex(Boundary_Maps, Matrix_Representation_Embbedings)
    Down = DownLaplacians_Of_Complex(Boundary_Maps, Matrix_Representation_Embbedings)
    Complex_LabLacians = [Up[i] + Down[i] for i in range(len(Boundary_Maps))]
    return Complex_LabLacians



def Adjoint_Operator(LinMap, Embedding_Source, Embedding_Target):
    """
      Compute the adjoint operator of a linear map.

      Args:
          LinMap (array): Linear map as an array.
          Embedding_Source (array): Embedding of the source space.
          Embedding_Target (array): Embedding of the target space.

      Returns:
          array: The adjoint operator of the linear map.
      """
    T = np.array(LinMap).astype('float64')
    Is = np.array(Embedding_Source).astype('float64')
    It = np.array(Embedding_Target).astype('float64')
    # E = B.conjugate().T # adjoint operator
    if np.abs(np.linalg.det(Is.T@Is)) < 0.001: # if the source space is zero space
      E = np.zeros_like(T).T
    else:
      Qs = np.linalg.inv((Is.conjugate().T@Is))
      Qt = (It.conjugate().T@It)
      E = Qs@(T.conjugate().T)@Qt # adjoint operator
    return E



###################################################
###################################################



def Persistent_Laplacian(Complex_Basis_Reps_Matrices, Persistence_Complex_Reps_Matrices, PCn, Tn):
    """
      Compute the persistent Laplacian matrices.

      Args:
          Complex_Basis_Reps_Matrices (list): List of complex basis representation matrices.
          Persistence_Complex_Reps_Matrices (list): List of persistence complex representation matrices.
          PCn (list): List of boundary maps.
          Tn (list): List of maps.

      Returns:
          list: List of persistent Laplacian matrices.
      """
    AReps = Complex_Basis_Reps_Matrices
    CReps = Persistence_Complex_Reps_Matrices
    L = len(PCn)
    Ln = [[] for i in range(L)]

    for i in range(1, L-1):
      B = np.array(PCn[i]).astype(np.float64)
      B_hat = Adjoint_Operator(B, AReps[i], AReps[i-1])
      DownLap = B_hat@B
      G = np.array(Tn[i+1]).astype(np.float64)
      G_hat = Adjoint_Operator(G, CReps[i+1], AReps[i])
      UpLap = G@G_hat
      Ln[i] = DownLap + UpLap

    G = np.array(Tn[1]).astype(np.float64)
    G_hat = Adjoint_Operator(G, CReps[1], CReps[0])
    UpLap = G@G_hat
    Ln[0] = UpLap

    B = np.array(PCn[L-1]).astype(np.float64)
    B_hat = Adjoint_Operator(B, AReps[L-1], AReps[L-2]) # in our cases, L>1, or at least make the complex not just A_0, otherwise ridiculous
    DownLap = B_hat@B
    Ln[L - 1] = DownLap

    return Ln



def ComplexFeatures(Complex_Basis_Reps_Matrices, Complex_Boundary_Matrices, Length):
    """
      Compute features of a complex.

      Args:
          Complex_Basis_Reps_Matrices (list): List of complex basis representation matrices.
          Complex_Boundary_Matrices (list): List of complex boundary matrices.
          Length (int): Length of the complex.

      Returns:
          list: List of complex features.
      """
    Ups = UpLaplacians_Of_Complex(Complex_Boundary_Matrices, Complex_Basis_Reps_Matrices)
    Downs = DownLaplacians_Of_Complex(Complex_Boundary_Matrices, Complex_Basis_Reps_Matrices)
    Laps = Laplacians_Of_Complex(Complex_Boundary_Matrices, Complex_Basis_Reps_Matrices)

    True_Zero_Space = [1 for i in range(Length)]
    for i in range(Length):
      if sympy.det(sympy.Matrix(Complex_Basis_Reps_Matrices[i]).T*sympy.Matrix(Complex_Basis_Reps_Matrices[i])) < 0.00001:
        True_Zero_Space[i] = 0

    Complex_Features = {}
    Complex_Features.update({'Ups_Features': [MatrixFeatures(Ups[i]) for i in range(len(Ups))]})
    Complex_Features.update({'Downs_Features': [MatrixFeatures(Downs[i]) for i in range(len(Downs))]})
    Complex_Features.update({'Laps_Features': [MatrixFeatures(Laps[i]) for i in range(len(Laps))]})

    for i in range(Length):
      Complex_Features['Ups_Features'][i]['Null'] = True_Zero_Space[i]*Complex_Features['Ups_Features'][i]['Null']
      Complex_Features['Downs_Features'][i]['Null'] = True_Zero_Space[i]*Complex_Features['Downs_Features'][i]['Null']
      Complex_Features['Laps_Features'][i]['Null'] = True_Zero_Space[i]*Complex_Features['Laps_Features'][i]['Null']

    DownsNulls = [Complex_Features['Downs_Features'][i]['Null'] for i in range(len(Complex_Features['Downs_Features']))]
    LapsNull = [Complex_Features['Laps_Features'][i]['Null'] for i in range(len(Complex_Features['Laps_Features']))]
    DownsEigens = [Complex_Features['Downs_Features'][i]['Eigens'] for i in range(len(Complex_Features['Downs_Features']))]
    LapsEigens = [Complex_Features['Laps_Features'][i]['Eigens'] for i in range(len(Complex_Features['Laps_Features']))]

    Complex_Features.update({'Dirac_Features': Dirac_Features(DownsNulls, LapsNull, DownsEigens, LapsEigens)})
    # Complex_Features.update({'Dirac_Nullity': Dirac_Nullity(LapsNull, DownsNulls)})


    return [Complex_Features, Ups, Downs, Laps]

def Persistence_Features(Complex_Basis_Reps_Matrices1, Persistence_Complex_Reps_Matrices, Complex_Boundary_Matrices1, Complex_Diagonal_Matrices2):
    """
      Compute features of persistence complexes.

      Args:
          Complex_Basis_Reps_Matrices1 (list): List of complex basis representation matrices.
          Persistence_Complex_Reps_Matrices (list): List of persistence complex representation matrices.
          Complex_Boundary_Matrices1 (list): List of complex boundary matrices.
          Complex_Diagonal_Matrices2 (list): List of complex diagonal matrices.

      Returns:
          list: List of persistence features.
      """
    Pers_Laps = Persistent_Laplacian(Complex_Basis_Reps_Matrices1, Persistence_Complex_Reps_Matrices, Complex_Boundary_Matrices1, Complex_Diagonal_Matrices2)
    Pers_Features = [MatrixFeatures(Pers_Laps[i]) for i in range(len(Pers_Laps))]
    True_Zero_Space = [1 for i in range(len(Pers_Laps))]

    for i in range(len(Pers_Laps)):
      if sympy.det(sympy.Matrix(Complex_Basis_Reps_Matrices1[i]).T*sympy.Matrix(Complex_Basis_Reps_Matrices1[i])) < 0.00001:
        True_Zero_Space[i] = 0

    for i in range(len(Pers_Laps)):
      Pers_Features[i]['Null'] = True_Zero_Space[i]*Pers_Features[i]['Null']
    return [Pers_Features, Pers_Laps]





