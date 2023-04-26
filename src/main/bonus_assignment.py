import numpy as np

np.set_printoptions(precision=7, suppress=True, linewidth=100)

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

def absolute_error(precise, approximate): 
    op = precise-approximate
    return abs(op)

def gauss_seidel_update(A, b, x, tolerance):

    stop = False
    count = 0

    while (stop != True):
        
        for n in range(len(A)):
            subs = 0

            for m in range(len(A)):
                if (n!=m):
                    subs -= A[n][m] * x[m]

            new = (b[n] - subs ) / A[n][n]
            if (absolute_error(x[n], new) <= tolerance ):
                stop = True
                break
            x[n] = new
            
        count = count + 1
 
    return count

# xn = (bn - a1nx1 - a2nx2 - ... - a(n-1)n x(n-1)) / ann
def jacobi_update(A, b, x, tolerance):

    stop = False
    count = 0

    while (stop != True):
        x_copy = np.copy(x)

        for n in range(len(A)):
            subs = 0

            for m in range(len(A)):
                if (n!=m):
                    subs -= A[n][m] * x_copy[m]

            new = (b[n] - subs ) / A[n][n]
            
            if (absolute_error(x[n], new) <= tolerance ):
                stop = True
                break

            x[n] = new
            
        count = count + 1
        
    return count


def function(t: float, w: float):
    return w - t*t*t


def do_work(t, w, h):
    basic_function_call = function(t, w)

    incremented_t = t + h
    incremented_w = w + (h * basic_function_call)
    incremented_function_call = function(incremented_t, incremented_w)

    return basic_function_call + incremented_function_call

def modified_eulers():
    original_w = .5
    start_of_t, end_of_t = (0, 3)
    num_of_iterations = 100

    # set up h
    h = (end_of_t - start_of_t) / num_of_iterations

    for cur_iteration in range(0, num_of_iterations):
        # do we have all values ready?
        t = start_of_t
        w = original_w
        h = h

        # create a function for the inner work
        inner_math = do_work(t, w, h)

        # this gets the next approximation
        next_w = w + ( (h / 2) * inner_math )

        # we need to set the just solved "w" to be the original w
        # and not only that, we need to change t as well
        start_of_t = t + h
        original_w = next_w
        
    return original_w

def custom_derivative(value):
    return (3 * value* value) - (2 * value)

def newton_raphson(initial_approximation: float, tolerance: float, sequence: str):
    # remember this is an iteration based approach...
    iteration_counter = 0

    # finds f
    x = initial_approximation
    f = eval(sequence)

    # finds f' 
    f_prime = custom_derivative(initial_approximation)
    
    approximation: float = f / f_prime
    while(abs(approximation) >= tolerance):
        # finds f
        x = initial_approximation
        f = eval(sequence)

        # finds f' 
        f_prime = custom_derivative(initial_approximation)

        # division operation
        approximation = f / f_prime

        # subtraction property
        initial_approximation -= approximation
        iteration_counter += 1

    return iteration_counter

def apply_div_diff(matrix):

    size = len(matrix)
    for i in range(2, size):
        for j in range(2,i+2):
            
            if j >= len(matrix[i]) or matrix[i][j] != 0:
                continue

            # something get left and diag left
            left = matrix[i][j-1]
            diag_left = matrix[i-1][j-1]
            numerator = left - diag_left

            denominator = matrix[i][0] - matrix[i-j+1][0]

            operation = numerator / denominator
            matrix[i][j] = operation
    return matrix


def hermite_interpolation(x_points, y_points, slopes):

    #main difference with hermite's method , using instances with x 

    num_of_points = len(x_points)
    matrix = np.zeros((num_of_points * 2, num_of_points * 2))

    #populate x values

    index = 0
    for x in range (0, len(matrix), 2):
        matrix[x][0]= x_points[index]
        matrix[x+1][0] = x_points[index]
        index += 1

    # prepopulate y values
    index = 0
    for x in range (0, len(matrix), 2):
        matrix[x][1]= y_points[index]
        matrix[x+1][1] = y_points[index]
        index += 1

    #prepopulate with derivatives (every other row)
    index = 0
    for x in range(1, len(matrix), 2):
        matrix[x][2] = slopes[index]
        index += 1

    apply_div_diff(matrix)
    print(matrix)

# jacobi / gauss-steidel
# define the equations ahead of time
    # know your pivots (a_ii)
    # know what other coefficients are in that row
    # pass in the x _values
# save x_n and repeat for a given threshold

if __name__ == "__main__":

    init_guess = np.array([0, 0,0], dtype=np.double)
    tolerance = 0.000001
    matrix = np.array([[3, 1, 1],
                       [1, 4, 1],
                       [2, 3, 7]], dtype=np.double)
    b_vector = np.array([1, 3, 0], dtype=np.double)


    newA1, newb1 = make_diagonally_dominant(matrix, b_vector)
    newA2, newb2 = make_diagonally_dominant(matrix, b_vector)

    print(gauss_seidel_update(newA1, newb1, init_guess, tolerance), end="\n\n")

    print(jacobi_update(newA2, newb2, init_guess, tolerance), end="\n\n")
    print(newton_raphson(0.5, tolerance, "x*x*x - x*x + 2"), end="\n\n")

    x_points = [0, 1, 2]
    y_points =  [1, 2, 4]
    slopes = [1.06, 1.23, 1.55]
    hermite_interpolation(x_points, y_points, slopes)
    print("\n")

    print(round(modified_eulers(),5))