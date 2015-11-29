export view_specs, bounding_box

import Requests
import FileIO, Images, Colors

const REL_DATA_DIR = joinpath("..", "data")
if !isdir(REL_DATA_DIR)
  mkdir(REL_DATA_DIR)
end
const DATA_DIR = Base.realpath(REL_DATA_DIR)
const BASE_URL = "http://vision.middlebury.edu/mview/data/data/"

const DATA_SETS = Dict{Symbol,Tuple{AbstractString,AbstractString,Bool}}(
    :dino               => ("dino"            , "dino"    , false),
    :dino_ring          => ("dinoRing"        , "dinoR"   , false),
    :dino_sparse_ring   => ("dinoSparseRing"  , "dinoSR"  , false),
    :temple             => ("temple"          , "temple"  , true ),
    :temple_ring        => ("templeRing"      , "templeR" , true ),
    :temple_sparse_ring => ("templeSparseRing", "templeSR", true )
)

function check_dataset(ds_key)
  (fn_base, im_fn_base, is_rgb) = DATA_SETS[ds_key]

  ds_dir = joinpath(DATA_DIR, fn_base)
  if !isdir(ds_dir)
    info("$ds_dir not found.")
    zip_file = joinpath(DATA_DIR, fn_base*".zip")
    zip_url = BASE_URL*fn_base*".zip"

    if !isfile(zip_file)
      info("$zip_file not found. Fetching from $zip_url...")
      Requests.save(Requests.get(zip_url), zip_file)
    end

    run(`unzip $zip_file -d $DATA_DIR`)
  end
end

function view_specs(ds)
  check_dataset(ds)
  data_set = DATA_SETS[ds]
  dir = joinpath(DATA_DIR, data_set[1])
  local lines
  open(joinpath(dir, data_set[2] * "_par.txt")) do f
    readline(f)
    lines = readlines(f)
  end

  function ln_to_spec(ln)
    strs = split(chomp(ln))
    name = strs[1]
    nums = map(s -> parse(Float32, s), strs[2:end])
    K = reshape(nums[1:9],(3,3))'
    R = reshape(nums[10:18],(3,3))'
    T = nums[19:21]
    cam = M34(K * hcat(R,T))
    (joinpath(dir, name), cam, data_set[3])
  end

  map(ln_to_spec, lines)
end

function bounding_box(ds)
  check_dataset(ds)

  function line_to_vec(ln)
    V3(map(float, split(chomp(ln)[2:end-1])))
  end

  data_set = DATA_SETS[ds]
  local lo, hi
  open(joinpath(DATA_DIR, data_set[1], "README.txt")) do f
    while !eof(f)
      ln = readline(f)
      if ismatch(r"^The \(tight\) bounding box", ln)
        lo = line_to_vec(readline(f))
        hi = line_to_vec(readline(f))
        break
      end
    end
  end
  vcat(lo',hi')
end

function View(spec)

  # read in image
  image = FileIO.load(spec[1])
  (nx, ny) = size(image)

  # convert to float array with alpha channel
  local im
  if spec[3]
    im = Array(Float32, (4, nx, ny))
    for i in 1:nx, j in 1:ny
      im[1,i,j] = 1.0
      samp = image[i,j]
      im[2,i,j] = Float32(samp.r)
      im[3,i,j] = Float32(samp.g)
      im[4,i,j] = Float32(samp.b)
    end
  else
    image = convert(Images.Image{Colors.Gray}, image)
    im = Array(Float32, (2, nx, ny))
    for i in 1:nx, j in 1:ny
      im[1,i,j] = 1.0
      im[2,i,j] = Float32(image[i,j])
    end
  end

  View(spec[2], im)
end
