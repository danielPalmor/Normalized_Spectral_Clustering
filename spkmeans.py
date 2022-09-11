import sys

import numpy as np
import pandas as pd
import spectralkmeans


def normalized_spectral_clustering(k, goal, file_name):
    """
    This function performs the full normalized spectral clustering algorithm
    :param k: The numbers of clusters
    :param goal: The goal
    :param file_name: The path to the file that contains the points
    :return: 0 for success, 1 for failure
    """
    np.random.seed(0)
    X = pd.read_csv(file_name, header=None)
    rows, columns = np.shape(X.values)
    X = np.reshape(X.values, newshape=X.size)

    X_as_list = [scalar for scalar in X]

    if goal == "spk":

        result = spectralkmeans.fit(X_as_list, 1, rows, columns, k)
        k, result = int(result[-1]), result[:-1]
        T_as_df = pd.DataFrame(np.array(result).reshape(rows, k))

        if T_as_df.empty or k >= len(T_as_df):
            print("An Error Has Occurred")
            return 1

        # Performing kmeans++
        rand_indexes = pd.Series(T_as_df.index.values)  # The indexes to choose from randomly
        index = np.random.choice(T_as_df.index.values)  # The chosen index
        indexes = [rand_indexes[index]]  # The list of chosen indexes

        first_centroid = T_as_df.iloc[index]
        initial_centroids = [np.array(first_centroid)]
        i = 1
        while i < k:
            list_of_p = np.array([])
            for point in T_as_df.values:
                # Calculating the min distance of the current point from all centroids
                d_point = min(np.dot(point - centroid, point - centroid)
                              for centroid in initial_centroids)
                list_of_p = np.append(list_of_p, d_point)
            list_of_p = list_of_p / np.sum(list_of_p)  # Calculating the probabilities
            i += 1

            index = np.random.choice(T_as_df.index.values, p=list_of_p)
            indexes.append(rand_indexes[index])
            next_centroid = T_as_df.iloc[index]
            initial_centroids.append(next_centroid.values)

        # Adjustments for the CAPI
        rows, columns = np.shape(initial_centroids)
        initial_centroids_as_list = np.reshape(initial_centroids, newshape=rows * columns)
        initial_centroids_as_list = [scalar for scalar in initial_centroids_as_list]

        the_centroids = spectralkmeans.kmeans_fit(k, result,
                                                  initial_centroids_as_list,
                                                  len(initial_centroids[0]),
                                                  len(T_as_df.values))
        # Printing the chosen indexes
        for i in range(len(indexes)):
            if i == len(indexes) - 1:
                print(int(indexes[i]))
            else:
                print(int(indexes[i]), end=",")

        # Printing the chosen centroids
        for i, point in enumerate(the_centroids):
            if (i + 1) % columns == 0:
                print(format(point, ".4f"))
            else:
                print(format(point, ".4f"), end=",")
        return 0
    elif goal == "wam":
        result = spectralkmeans.fit(X_as_list, 2, rows, columns, k)
    elif goal == "ddg":
        result = spectralkmeans.fit(X_as_list, 3, rows, columns, k)
    elif goal == "lnorm":
        result = spectralkmeans.fit(X_as_list, 4, rows, columns, k)
    elif goal == "jacobi":
        result = spectralkmeans.fit(X_as_list, 5, rows, columns, k)
    else:
        print("Invalid Input!")
        return 1

    for i, scalar in enumerate(result):
        if (i + 1) % rows == 0:
            print(format(scalar, ".4f"))
        else:
            print(format(scalar, ".4f"), end=",")
    return 0


def main(argv):
    if len(argv) != 4 or not argv[1].isnumeric():
        print("Invalid Input!")
        return 1

    if int(argv[1]) < 0:
        print("Invalid Input!")
        return 1

    normalized_spectral_clustering(int(argv[1]), argv[2], argv[3])
    return 0


if __name__ == "__main__":
    main(sys.argv[:])
