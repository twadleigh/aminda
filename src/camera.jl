export View, homogenize, project, homography, bounds, camera_resolution

type View
  camera :: M34
  image :: Array{Float32, 3}
end

homogenize(v :: V3) = V4(v[1], v[2], v[3], 1.0)
*(m :: M4, v :: V3) = m * homogenize(v)
*(m :: M34, v :: V3) = m * homogenize(v)
project(v :: V4) = V3(v[1]/v[4], v[2]/v[4], v[3]/v[4])
project(v :: V3) = V2(v[1]/v[3], v[2]/v[3])

"""
Given 2 camera matrices, calculate the homography that takes pixels from camera 1 to
camera 2 through camera 1's z-depth image plane.
"""
homography(c1, c2, z) = M3((Array(c2)/[Array(c1) ; 0 0 0 z])*[eye(3) ; 0 0 1])

function corners(voi)
  (lo,hi) = voi

  [ V3(lo[1], lo[2], lo[3])
    V3(hi[1], lo[2], lo[3])
    V3(lo[1], hi[2], lo[3])
    V3(hi[1], hi[2], lo[3])
    V3(lo[1], lo[2], hi[3])
    V3(hi[1], lo[2], hi[3])
    V3(lo[1], hi[2], hi[3])
    V3(hi[1], hi[2], hi[3])
  ]
end

function bounds(cam, voi)

  function _proj(p)
    px = cam * p
    V3(px[1]/px[3], px[2]/px[3], px[3])
  end

  crnrs = corners(voi)
  xyzs = zip(map(_proj, crnrs)...)
  mns = map(minimum, xyzs)
  mxs = map(maximum, xyzs)

  (V3(mns), V3(mxs))
end

# smallest resolvable distance at unit depth (roughly, 1/f)
camera_resolution(cam) = minimum(svdvals([Array(cam) ; 0 0 0 1]\eye(4,2)))
