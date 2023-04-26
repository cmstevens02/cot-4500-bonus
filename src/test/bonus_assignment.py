import numpy as np


# taken from Prof Parra's Notes
def make_diagonally_dominant(matrix, b_vector):
    n = len(matrix)

    for i in range(n):
        pivot: float = matrix[i][i]
        sum_of_other_elements = sum(abs(matrix[i][i+1:]))

        # we can guarantee this pivot is the largest in the row
        if abs(pivot) > abs(sum_of_other_elements):
            continue

        # if we reach this point, this means we need to swap AT LEAST ONCE
        max_value_of_row = 0
        max_index_in_row = 0
        for j in range(n):
            current_value_in_row: float = abs(matrix[i][j])
            if current_value_in_row > max_value_of_row:
                max_value_of_row = current_value_in_row
                max_index_in_row = j

        # now that we have a new "pivot", we swap cur_row with the expected index
        matrix[[i, max_index_in_row]] = matrix[[max_index_in_row, i]]
        b_vector[[i, max_index_in_row]] = b_vector[[max_index_in_row, i]]
        
    return matrix, b_vector

# xn = (bn - a1nx1 - a2nx2 - ... - a(n-1)n x(n-1)) / ann
def jacobi_update(A, b, x):

    for n in range(len(A)):
        subs = 0; 

        for m in range(n):
            subs -= a[n][m] * m 
        
        x[n] = b[n] - subs 


def jacobi(init, tol, A, b):
    
    count = 0

    while (diagDomCheck(A, b) == False):
        jacobi_update(A, b, init)
        count += 1

    print(count)



def diagDomCheck(matrix, b):
    A = np.concatenate((matrix, b.reshape(n, 1)), axis=1)

    length = len(A)

    # the question just say diagonally dominant, not strictly so...
    for i in range(length):
        diag = A[i][i]
        others = np.sum(A[i]) - diag

        if (diag < others):
            return False

    return True
    

# jacobi / gauss-steidel
# define the equations ahead of time
    # know your pivots (a_ii)
    # know what other coefficients are in that row
    # pass in the x _values
# save x_n and repeat for a given threshold

if __name__ == "__main__":

    init_guess = np.array([0, 0,0])
    tolerance = 0.000001
    matrix = np.array([[3, 1, 1],
                       [1, 4, 1],
                       [2, 3, 7]])
    b_vector = np.array([1, 3, 0])

    jacobi(init_guess, tolerance, matrix, b_vector)

