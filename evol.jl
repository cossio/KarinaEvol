import ArgParse, Distributions

const arg_parse_settings = ArgParse.ArgParseSettings()

@ArgParse.add_arg_table arg_parse_settings begin
	"--g"
		help = "rate de crecimiento clon inicial del tumor"
		arg_type = Float64
		range_tester = g -> g ≥ 0
		required = true
	"--d"
		help = "rate muerte clon inicial del tumor"
		arg_type = Float64
		range_tester = d -> d ≥ 0
		required = true
	"--kill"
		help = "kill clon inicial del tumor"
		arg_type = Float64
		range_tester = kill -> kill ≥ 0
		required = true
	"--lambda"
		help = "lambda clon inicial del tumor"
		arg_type = Float64
		range_tester = λ -> λ ≥ 0
		required = true
	"--si"
		help = "condiciones iniciales SI"
		arg_type = String
		range_tester = isfile
		required = true
	"--mutfreq"
		help = "frecuencia de mutaciones por division"
		arg_type = Float64
		range_tester = f -> 0 ≤ f < Inf
		required = true
	"--out"
		help = "fichero salida"
		arg_type = String
		required = true
	"--T"
		help = "tiempo de corrida"
		arg_type = Float64
		default = 1000.
		range_tester = T -> 0 <= T <= Inf
end

const args = ArgParse.parse_args(arg_parse_settings)

const γ = 8.0
const Npeptidos = 20

const s = 5		# sitios APC
const ϵ = 21	# rate proliferacion Eff (CHEQUEAR este valor)
const μ = 9		# rate proliferacion Reg (CHEQUEAR)
const Am = 1e4	# maximo de APC expandidos por tumor

const h = 2 	# exponencial en calculo afinidad total del clon T por APC
const Km = 16	# maxima afinidad de T cel por APC  (CHEQUEAR)

# limits for conjugation constants in thymic selection (Kalet2003)
const AFmi = 0.08
const AFmx = 8

# limits for random generation of affinity
const kmi = (AFmx / (Km - AFmx))^(1/h) / 0.7
const kmx = (AFmi / (Km - AFmi))^(1/h) / 0.001
@assert kmi < kmx
# valores en Kalet2003 (intercambio TPmi con TPmx para q el menor valor quede como TPmi=1.42 y el mayor como TPmx=55.9)

const T = 1000.0

const Sc = 24 / log(2)	# Tcell clones generated per unit time
const renovate_Δt = 1e-1 / Sc	# intervalo d tiempo entre renovates (tumor y SI)

#@assert false	# para asegurarte de chequear los numeritos arriba, princesa. Quitar esto para correr. Besos, tu shrek.


# a clone of tumor cells
type Tumor
	g::Float64
	λ::Float64
	d::Float64
	kill::Float64
	id::Int
	parent::Int
	Tu::Float64
	pep::Vector{Float64}

	# Runge-Kutta stuff
	Ktu::Vector{Float64}
	Xtu::Vector{Float64}

	function Tumor(g::Float64, λ::Float64, d::Float64, kill::Float64, id::Int, parent::Int, Tu::Float64, pep::Vector{Float64})
		@assert g ≥ 0
		@assert λ ≥ 0
		@assert d ≥ 0
		@assert kill ≥ 0
		@assert Tu ≥ 0
		@assert all(pep .≥ 0)
		@assert length(pep) == Npeptidos
		new(g, λ, d, kill, id, parent, Tu, pep, zeros(6), zeros(6))
	end
end

# a clone of T cells (eff and reg)
type TCell
	e::Float64
	r::Float64
	k::Vector{Float64}	# aff por peptidos
	Aff::Float64	# aff por APC

	# Runge-Kutta stuff
	Ke::Vector{Float64}
	Kr::Vector{Float64}
	Xe::Vector{Float64}
	Xr::Vector{Float64}

	function TCell(e::Float64, r::Float64, k::Vector{Float64}, Aff::Float64)
		@assert e ≥ 0
		@assert r ≥ 0
		@assert all(k .≥ 0)
		@assert length(k) == Npeptidos
		@assert Aff ≥ 0
		new(e, r, k, Aff, zeros(6), zeros(6), zeros(6), zeros(6))
	end
end

# T cells at tumor site of a clone
et(ec::Float64, Rc::Float64, a::Float64) = ϵ * (1 - Rc / (s * a))^(s - 1) * ec
rt(rc::Float64, Ec::Float64, a::Float64) = μ * (s - 1) / (s * a) * rc * Ec

function APC_a(tumor::Set{Tumor})
	imm = isempty(tumor) ? 0. : sum(t.λ * t.Tu for t in tumor)
	return 1 + Am * imm / (Am + imm)
end

# T cells bound at APC. Where f is APC free sites
frac_bound(si::TCell, free::Float64, vol::Float64) = si.Aff * free / (vol + si.Aff * free)
ec(si::TCell, free::Float64, vol::Float64) = frac_bound(si, free, vol) * si.e
rc(si::TCell, free::Float64, vol::Float64) = frac_bound(si, free, vol) * si.r

function Aff!(SI::Set{TCell}, tumor::Set{Tumor}, self::Vector{Float64})
	if !isempty(tumor); Ptot = 1. + sum(t.λ * t.Tu * sum(t.pep) for t in tumor); end
	# 1 = sum(self), pq ya esta normalizado
	for si in SI
		σ = self ⋅ si.k
		if !isempty(tumor)
			σ += sum(t.λ * t.Tu * t.pep ⋅ si.k for t in tumor)
			σ /= Ptot
		end
		si.Aff = Km * σ^h / (1 + σ^h)
	end
end

# APC free sites
function free(SI::Set{TCell}, tumor::Set{Tumor}, a::Float64)
	vol = sum(si.e + si.r for si in SI) + a
	# pto fijo de Newton para hallar f (total free sites en APC adimensionalizado)
	# resolver h(f) = 0
	h(f)  = sum(ec(si, f, vol) + rc(si, f, vol) for si in SI) + f - s * a
	dh(f) = sum(si.Aff * (si.e + si.r) * vol / (vol + si.Aff * f)^2 for si in SI) + 1
	Δ = Inf
	f = a * s / 2
	while Δ > 1e-6
		fnew = f - h(f) / dh(f)
		Δ = abs(fnew - f)
		f = fnew
	end
	@assert f ≥ 0
	return f
end

# genera afinidades aleatorias para celulas T
function rand_k!(k::Vector{Float64})
	@assert length(k) == Npeptidos == 20
	for l = 1:10
		if l ≤ 10 && rand() < 0.5 || l > 10 && rand() < exp10((i - 10) / 5 - 2)
			k[l] = rand(Distributions.Uniform(kmi, kmx))
		else
			k[l] = 0
		end
	end
end

function randTCell(self::Vector{Float64})
	k = zeros(Npeptidos)
	Aff = 0.0
	while !(AFmi ≤ Aff ≤ AFmx)
		rand_k!(k)
		σ = self ⋅ k
		Aff = Km * σ^h / (1 + σ^h)	# Aff por self solo!
	end
	return TCell(0.0017, 0.041, k, Aff)
end

function renovate!(SI::Set{TCell}, self::Vector{Float64})
	# quita clones T con pocas celulas
	filter!(si -> max(si.e, si.r) ≥ 0.02, SI)
	# genera clones nuevos
	clones_nuevos = rand(Distributions.Poisson(Sc * renovate_Δt))
	for i = 1 : rand(Distributions.Poisson(Sc * renovate_Δt))
		push!(SI, randTCell(self))
	end
end

function randpar(μ::Float64, min::Float64, max::Float64, cv::Float64 = 5.0)
	@assert μ > 0 && max > 0 && min > 0 && min < max && cv > 0
	exp10(rand(Distributions.TruncatedNormal(log10(μ), (log10(max) - log10(min)) / cv, log10(min), log10(max))))
end

function rand_clon(parent::Tumor, id::Int)
	g = randpar(parent.g, 0.01, 42.)
	λ = randpar(parent.λ, 1e-5, 1000.)
	d = parent.d > 0 ? randpar(parent.d, 0.01, 42.) : parent.d
	kill = randpar(parent.kill, 4.2, 420.)
	return Tumor(g, λ, d, kill, id, parent.id, 1e-4, rand_pep(parent.pep))
end

function renovate!(tumor::Set{Tumor}, last_id::Int)
	# quita clones tumorales con pocas celulas
	filter!(clon -> clon.Tu ≥ 0.0001, tumor)

	# generar clones nuevos mutados
	childs = Set{Tumor}()
	for clon in tumor
		# probability of mutation
		mutprob = clon.g * clon.Tu * args["mutfreq"] * renovate_Δt
		for k = 1 : rand(Distributions.Poisson(mutprob))
			last_id += 1
			newclone = rand_clon(clon, last_id)
			push!(childs, newclone)
		end
	end
	isempty(childs) || push!(tumor, childs...)
	return last_id
end

function rand_pep(parent_pep::Vector{Float64})
	pep = randpar.([parent_pep[1:10]; fill(1., 10)], 1e-4, 500.)
	pep ./= sum(pep)
	return pep
end

# Runge–Kutta–Fehlberg method of order 4
function RGK_5_6_step!(SI::Set{TCell}, tumor::Set{Tumor}, self::Vector{Float64}, Δt::Float64)

	# esta funcion solo calcula los K de los 6 stages del Runge-Kutta

	const tab::Matrix{Float64} = [0			0			0			0			0		0;	# coeffs k1
				 				  1/4		0			0			0			0		0;	# coeffs k2
	             				  3/32		9/32		0			0			0		0;	# coeffs k3
				 			  	  1932/2197	-7200/2197	7296/2187	0			0		0;	# coeffs k4
				 			  	  439/216	-8			3680/513	-845/4104	0		0;	# coeffs k5
				 			  	  -8/27		2			-3544/2565	1859/4104	-11/40	0;]	# coeffs k6

	# six stages of Runge Kutta
	for i = 1:6
		for clon in tumor
			if i == 1
				clon.Xtu[i] = clon.Tu
			else
				clon.Xtu[i] = clon.Tu = clon.Xtu[1] + sum(tab[i,j] * clon.Ktu[j] for j = 1 : i - 1)
			end
		end

		for si in SI
			if i == 1
				si.Xe[i] = si.e
				si.Xr[i] = si.r
			else
				si.Xe[i] = si.e = si.Xe[1] + sum(tab[i,j] * si.Ke[j] for j = 1 : i - 1)
				si.Xr[i] = si.r = si.Xr[1] + sum(tab[i,j] * si.Kr[j] for j = 1 : i - 1)
			end
		end

		# updatea afinidades TCell/APC
		Aff!(SI, tumor, self)
		# updatea otros parametros
		a = APC_a(tumor)
		f = free(SI, tumor, a)	# free conjugation sites
		volLy = sum(si.e + si.r for si in SI) + a		# vol Lymph compartment

		Ec = isempty(SI) ? 0. : sum(ec(si,f, volLy) for si in SI)	# total conjugated E cells
		Rc = isempty(SI) ? 0. : sum(rc(si,f, volLy) for si in SI)	# total conjugated R cells

		Et = isempty(SI) ? 0. : sum(et(ec(si,f,volLy), Rc, a) for si in SI)	# total E cells at tumor site
		Rt = isempty(SI) ? 0. : sum(rt(rc(si,f,volLy), Ec, a) for si in SI)	# total R cells at tumor site

		Tu = isempty(tumor) ? 0. : sum(clon.Tu for clon in tumor)	# total tumor cells
		volTu = Tu + Et + Rt	# tumor compartment volume

		# calcula los K del Runge-Kutta
		for clon in tumor
			clon.Ktu[i] = Δt * (clon.g - clon.kill * Et / volTu * (1 - Rt / volTu)^(γ - 1)) * clon.Tu
		end

		for si in SI
			cloneEc = ec(si, f, volLy)
			cloneRc = rc(si, f, volLy)

			si.Ke[i] = Δt * (et(cloneEc, Rc, a) - si.e + cloneEc)
			si.Kr[i] = Δt * (rt(cloneRc, Ec, a) - si.r + cloneRc)
		end
	end
end

# Runge–Kutta–Fehlberg of order 4
function RGK_5_6!(SI::Set{TCell}, tumor::Set{Tumor}, self::Vector{Float64}, tf::Float64)
	const TOL = 1e-6
	const tab_xf    = [25/216	0	1408/2565	2197/4104	-1/5	0]
	const tab_error = [1/360	0	-128/4275	-2197/75240	1/50	2/55]

	Δt = renovate_Δt / 20		# step de integracion inicial
	t = 0.0

	while t < tf
		RGK_5_6_step!(SI, tumor, self, Δt)

		# error estimate
		ϵ = 0.0
		if !isempty(tumor)
			ϵ += sum((tab_error ⋅ clon.Ktu)^2 for clon in tumor) / Δt
		end
		if !isempty(SI)
			ϵ += sum((tab_error ⋅ si.Ke)^2 + (tab_error ⋅ si.Kr)^2 for si in SI) / Δt
		end

		if ϵ < TOL	# accept move?
			for clon in tumor
				clon.Tu += tab_xf ⋅ clon.Ktu
			end

			for si in SI
				si.e += tab_xf ⋅ si.Ke
				si.r += tab_xf ⋅ si.Kr
			end

			# updatea afinidades TCell/APC
			Aff!(SI, tumor, self)

			t += Δt
		end

		# step adaptativo
		δ = 0.84 * (TOL / ϵ)^0.25
		if δ < 0.1; Δt /= 10
		elseif δ > 4; Δt *= 4
		elseif 0.1 < δ < 4; Δt *= δ; end
		Δt = min(Δt, renovate_Δt / 2, (1 + 1e-10) * (tf - t))
	end
end


function main()
	self::Vector{Float64} = [10.0; 10.0; 0.001; fill(10.0, 6); 240.0; fill(0.001, 10)]
	self ./= sum(self)
	@assert length(self) == Npeptidos

	# inicializa T cells
	SI = Set{TCell}()

	SIinit = readdlm(args["si"])
	@assert size(SIinit, 2) == 2 + Npeptidos
	for i = 1 : size(SIinit, 1)
		e = SIinit[i, 1]
		r = SIinit[i, 2]
		k = SIinit[i, 2 + (1 : Npeptidos)]
		Aff = 0.
		push!(SI, TCell(e, r, k, Aff))
	end

	# inicializa tumor
	tumor = Set{Tumor}()
	last_tumor_id = 1
	push!(tumor, Tumor(args["g"], args["lambda"], args["d"], args["kill"], last_tumor_id, 0, 0.001, rand_pep(self)))


	out = open(args["out"], "w")

	# simular
	for t = 0.0 : renovate_Δt : args["T"]
		RGK_5_6!(SI, tumor, self, renovate_Δt)
		renovate!(SI, self)
		last_tumor_id = renovate!(tumor, last_tumor_id)

		# salva cosas
		r = isempty(SI) ? 0. : sum(si.r for si in SI)
		e = isempty(SI) ? 0. : sum(si.e for si in SI)
		tu = isempty(tumor) ? 0. : sum(clon.Tu for clon in tumor)
		println("$t\t$e\t$r\t$tu\t$(length(SI))\t$(length(tumor))")
		write(out, "$t\t$e\t$r\t$tu\t$(length(SI))\t$(length(tumor))\n")

		if tu > 1e5
			exit()
		end
	end

	close(out)

	# salvar estado equilibrio del SI.

	# eq = open("SI-eq.txt", "w")

	# for si in SI
	# 	write(eq, "$(si.e)\t$(si.r)")
	# 	for k in si.k
	# 		write(eq, "\t$k")
	# 	end
	# 	write(eq, "\n")
	# end

	# close(eq)

end


main()
