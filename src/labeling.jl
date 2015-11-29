export binary_labeling

# unregularized minimum cost binary labeling
function binary_labeling(σ)
  (nx,ny,nz) = size(σ)
  ϕ = similar(σ)
  for k in 1:nz, j in 1:ny, i in 1:nx
    ϕ[i,j,k] = σ[i,j,k] > 0 ? 0 : 1
  end
  ϕ
end

function binary_labeling(ρ, σ)
  const steps_per_test = 10
  const tol = 1.0e-4

  # initial guess is the unregularized estimate based on σ alone
  ϕ = binary_labeling(σ)

  # allocate workspace for iteration
  (nx,ny,nz) = size(σ)
  b = zeros(eltype(σ), (3,nx,ny,nz))
  d = similar(b)

  while true
    # update the labeling
    ϕₒ = copy(ϕ)
    LibAminda.update_labeling(steps_per_test, ρ, σ, ϕ, d, b)

    # test convergence
    if maximum(abs(ϕ-ϕₒ)) <= tol
      break
    end
  end

  ϕ
end
