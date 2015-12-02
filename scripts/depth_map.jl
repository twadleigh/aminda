push!(LOAD_PATH, realpath(pwd() * "/../src"))

using Aminda, ImageView

vs = view_specs(:dino_sparse_ring)
voi = bounding_box(:dino_sparse_ring)
vl = View(vs[1])
v = View(vs[2])
vr = View(vs[3])
nbrs = Set{View}([vl, vr])

dm = DepthMap(v, nbrs, voi)


view(dm.depth, pixelspacing = [1,1])
