module Aminda

module LibAminda
include("libaminda.jl")
end
import .LibAminda
include("types.jl")
include("camera.jl")
include("data.jl")
include("depth_map.jl")
include("labeling.jl")
include("fusion.jl")

end # module Aminda
