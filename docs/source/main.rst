main module
===========

.. automodule:: main
   :members:
   :undoc-members:
   :show-inheritance:

In this modules, the following steps were taken to perform PCA:

1. The data was read and converted to a np array.
2. The covariance of the genes was calculated using the :func:`covariance` function from the :mod:`AssignmentPCA` module.
3. From this matrix, the eigenvalues and eigenvectors were calculated and sorted by relevance.
4. The cumulative overall variance was then calculated.
5. A bar graph of principal component contribution and a scree plot were created.
6. It was decided to reduce to 3 dimensions.
7. A new subspace with only 3 PC's was created.
8. A 2D and 3D principal component plot were generated.
9. Loads for the three major PC's were calculated using :func:`calcLoads` from :mod:`AssignmentPCA`.
10. For each of the three major PC's, the 50 highest loads of the genes has been shown in a bar graph.
