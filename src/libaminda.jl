const LIBAMINDA = Base.source_dir() * "/libaminda.so"

function update_labeling(nits, ρ, σ, ϕ, d, b)
  (nx, ny, nz) = size(ϕ)
  ccall((:update_labeling, LIBAMINDA),
    Void, (Int32, Int32, Int32, Int32, Ptr{Float16}, Ptr{Float16}, Ptr{Float16}, Ptr{Float16}, Ptr{Float16}),
    nits, nx, ny, nz, ρ, σ, ϕ, d, b)
end

function map_to_plane(img, hom, nx, ny)
  @assert size(hom) == (3,3)
  (nc, inx, iny) = size(img)
  tgt = zeros(Float32, (nc, nx, ny))
  ccall((:map_to_plane, LIBAMINDA),
    Void, (Ptr{Float32}, Int32, Int32, Int32, Ptr{Float32}, Int32, Int32, Ptr{Float32}),
    hom', nc, inx, iny, img, nx, ny, tgt)
  return tgt
end

function mean_and_inverse_deviation(img, w)
  (nc, nx, ny) = size(img)
  mn = zeros(Float32, (nc, nx, ny))
  ccall((:mean_and_inverse_deviation, LIBAMINDA),
    Void, (Int32, Int32, Int32, Int32, Ptr{Float32}, Ptr{Float32}),
    w, nc, nx, ny, img, mn)
  return mn
end

function normalized_cross_correlation(im1, mn1, im2, mn2, w)
  (nc, nx, ny) = size(im1)
  nxcorr = fill(nx, ny, -1.0)
  ccall((:normalized_cross_correlation, LIBAMINDA),
    Void, (Int32, Int32, Int32, Int32, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}),
    w, nc, nx, ny, im1, mn1, im2, mn2, nxcorr)
  return nxcorr
end

function maximum(a,b)
  (nx, ny) = size(a)
  mx = Array(Float32, (nx, ny))
  ccall((:maximum, LIBAMINDA),
    Void, (Int32, Int32, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}),
    nx, ny, a, b, mx)
  return mx
end

function fill(nx, ny, x)
  a = Array(Float32, (nx, ny))
  ccall((:fill, LIBAMINDA),
    Void, (Int32, Int32, Ptr{Float32}, Float32),
    nx, ny, a, x)
  return a
end

function update_depth(nxc, mxnxc, zs, z)
  (nx, ny) = size(nxc)
  ccall((:update_depth, LIBAMINDA),
    Void, (Int32, Int32, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Float32),
    nx, ny, nxc, mxnxc, zs, z)
end
