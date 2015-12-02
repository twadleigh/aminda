export DepthMap

import ImageView

type DepthMap
  camera :: M34
  depth :: Array{Float32, 2}
  nxcorr :: Array{Float32, 2}
end

function DepthMap(view, nbrs, voi, w = 3)
  cam = view.camera
  im = view.image
  (nc, nx, ny) = size(im)
  mn = LibAminda.mean_and_inverse_deviation(im, w)

  # determine depth range, resolution
  bnds = bounds(cam, voi)
  near = bnds[1][3]
  far = bnds[2][3]
  dz = near * camera_resolution(cam)
  nz = ceil(Int, (far - near) / dz)
  dz = (far - near) / nz

  max_nxcorr = LibAminda.fill(nx, ny, -1.0)
  depth = LibAminda.fill(nx, ny, near)
  for k in 1:nz

    z = near + dz*(k-0.5)
    nxcorr = LibAminda.fill(nx, ny, -1.0)
    for nbr in nbrs
      cam2 = nbr.camera
      hom = Array(homography(cam, cam2, z))
      im2 = LibAminda.map_to_plane(nbr.image, hom, nx, ny)
      mn2 = LibAminda.mean_and_inverse_deviation(im2, w)
      nxc = LibAminda.normalized_cross_correlation(im, mn, im2, mn2, w)
      nxcorr = LibAminda.maximum(nxcorr, nxc)
    end

    LibAminda.update_depth(nxcorr, max_nxcorr, depth, z)
  end

  return DepthMap(cam, depth, max_nxcorr)
end
