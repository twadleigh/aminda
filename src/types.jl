using FixedSizeArrays

import Base.*

export M22, M23, M24, M32, M33, M34, M42, M43, M44, M2, M3, M4, V2, V3, V4

typealias M22 Mat{2, 2, Float32}
typealias M23 Mat{2, 3, Float32}
typealias M24 Mat{2, 4, Float32}
typealias M32 Mat{3, 2, Float32}
typealias M33 Mat{3, 3, Float32}
typealias M34 Mat{3, 4, Float32}
typealias M42 Mat{4, 2, Float32}
typealias M43 Mat{4, 3, Float32}
typealias M44 Mat{4, 4, Float32}

typealias M2 M22
typealias M3 M33
typealias M4 M44

typealias V2 Vec{2, Float32}
typealias V3 Vec{3, Float32}
typealias V4 Vec{4, Float32}
