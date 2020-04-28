def tetramesh(tri):
    import mpl_toolkits.mplot3d as a3
    import matplotlib.pyplot as plt
    axes = a3.Axes3D(plt.figure())
    # vts = xyz[tri, :]
    tri = a3.art3d.Poly3DCollection(dtri.points)
    tri.set_alpha(0.2)
    tri.set_color('grey')
    axes.add_collection3d(tri)
    axes.plot(point[:,0], point[:,1], point[:,2], 'ko')
    axes.set_axis_off()
    axes.set_aspect('equal')
    pl.show()