### A Pluto.jl notebook ###
# v0.19.13

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 7c2ceba1-6e8a-4a5f-845d-73b97810c099
begin
	using CSV, DataFrames, Query
	using StatsBase: mean, std
	using Optim, ForwardDiff
	using Turing, Distributions, MCMCChains
	using LinearAlgebra,PDMats
	using Plots, LaTeXStrings, Plots.Measures, StatsPlots
	using PlutoUI, PlutoTeachingTools
	using Downloads
	using Random
	using LombScargle
	using MarkdownLiteral
	Random.seed!(123)
end

# ╔═╡ f53f1d9f-8be8-4dcc-a637-f39b1eaacd63
md"""
**Astro 497 Project**
# Analyze Orbit Parameters
"""

# ╔═╡ 6d090107-d691-4a68-91ad-c7dae5df35ad
md"""
Our project aims to construct a dashboard that ingests radial velocity observation data and analyzes orbit parameters. In our dashboard, we
will utilize The California Legacy Survey (Rosenthal et al., 2021) data.
"""

# ╔═╡ 8cf8e22c-af0f-422b-89a1-109291cd749a
TableOfContents()

# ╔═╡ 8f829f16-553a-11ed-3fd9-9f77ddd265d5
md"""
# Ingest & select RV data to be analyzed
"""

# ╔═╡ a6eedc20-0171-48c7-938f-daee0ce4f8e9
begin
	fn = joinpath("../_assets/week4/legacy_data.csv")
	if !isfile(fn) || !(filesize(fn)>0)
		path = joinpath(pwd(),"data")
		mkpath(path)
		fn = joinpath(path,"legacy_data.csv")
		fn = Downloads.download("https://github.com/leerosenthalj/CLSI/raw/master/legacy_tables/legacy_data.csv", fn)
	end
	if filesize(fn) > 0
		df_all = CSV.read(fn, DataFrame)
		select!(df_all,["name","jd","mnvel","errvel", "cts","sval","tel"])
		# Rename columns to match labels from table in original paper
		rename!(df_all, "name"=>"Name","jd"=>"d","mnvel"=>"RVel","errvel"=>"e_RVel","tel"=>"Inst", "sval"=>"SVal")
		star_names = unique(df_all.Name)
		md"Successfully, read RVs for $(length(star_names)) stars from California Legacy Survey from [https://github.com/leerosenthalj/CLSI](https://github.com/leerosenthalj/CLSI) & from [Rosenthal et al. (2021)](https://doi.org/10.3847/1538-4365/abe23c) into `df_all`."
	else
		df_all = DataFrame()
		star_names = String[""]

		danger(md"Error reading data file with RVs.  Expect empty plots below.")
	end

end

# ╔═╡ e71bdd5e-a50e-46b2-8249-55ad7dffb789
begin
	select_star_id = PlutoRunner.currently_running_cell_id[] |> string
	select_star_url = "#$(select_star_id)"
	
md"""
Please select a star to be analyzed: $(@bind selected_star Select(collect(values(star_names)); default="117176"))
"""

end

# ╔═╡ 4a949edb-8737-4439-a1c0-14641bd99f8c
md"""
All available RV data for the selected star are plotted below.
"""

# ╔═╡ e4084b5c-19de-460d-aceb-ee6f3cfafbd9
md"""
Number of data points for each instrument is listed below. You can select the dataset to be analyzed based on this information.
"""

# ╔═╡ 64987161-afde-4b4d-93a8-19f7b40aeae7
md"""
The dashboard can handle data from one instrument or two instruments. The default option is to use two instruments. You can click the checkbox if you want to fit the models with data from only one instrument.\
Fit with only one instrument:$(@bind fit_with_one CheckBox())
"""

# ╔═╡ 0b49db68-973a-42e0-8eed-4f4405b5f808
md"""
# Fitting a linear model
"""

# ╔═╡ 5a929ea7-93cc-4c4d-9340-734bcd509719
md"""
A star without any planet is expected to have linear radial velocity. We fit a linear model to our data and see if it is a good fit. Sometimes, if the data seems to be linear, it is a good idea to subract linear part to see the deviations.
"""

# ╔═╡ 62e93800-22bd-40a1-9ee6-bcd1134539ae
md"""
The expected distribution of loss and the loss value at mle are plotted in the following figure.
"""

# ╔═╡ 09020d68-05d3-4908-9dcc-c1e0097cd54c
md"""
If you want to subtract linear fitting results from the data, click the checkbox: $(@bind subtract_linear CheckBox()). \
Warning: Only do this if the data looks to be linear and you want to investigate small deviations around the linear relationship. It is generally not suggested to subtract the linear fitting results.
"""

# ╔═╡ e53d6766-4fb1-4bda-822e-9e468b440fdf
md"""
# Lomb-Scargle
"""

# ╔═╡ 68a49dcb-3cd9-4905-abf6-d43083d256ce

warning_box(md"
1. Please try each instrument and use the one that gives the best periodogram.
1. Please visually inspect the peaks of power. To investigate different peaks (period), the user can change the periodogram range above. If there is a sharp peak with a very small period, it is suggested to set ~1% in the above box to detect possible hot Jupiters.
1. Periodiagram is not helpful for planets with periods longer than the observation duration.
")

# ╔═╡ 3b21b28a-52f0-4da0-beba-bb09d666b9c6
md"""
# Fitting one planet model
"""

# ╔═╡ 2384f27b-05d7-45f8-a18c-75bcdea14b80
md"""
The initial guesses are (P, K, h, k, M0-ω, C1, C2,σj):
"""

# ╔═╡ 6b5c4ec5-3b54-47e6-ad0b-d9a813b7bb7a
md"""
The results of model fitting are(P, K, h, k, M0-ω, C1, C2,σj):
"""

# ╔═╡ 69f83e46-8788-42ae-945d-3fb27b78a36f
md"""
To see how reasonable the fitted period is, check the Appendix for the plot where RV data is re-grouped by phase.
"""

# ╔═╡ 1d9f11ce-1515-4d75-9a30-55db3b3fa8ae
md"""
The RV observations used for the model fitting, together with the resulting model, are plotted below. 
"""

# ╔═╡ a494b98a-08df-464c-b3f4-584463f4210c
md"""
# Model Assessment
"""

# ╔═╡ fb7ddbe7-6f42-4a6d-8d28-611df2cdf549
@bind bootstrap confirm(PlutoUI.combine() do Child
md"""
**Bootstrap Parameter**

Number of samples  $(Child("num_bootstrap_samples", NumberField(10:10:10_000; default=100)))



Ready to run bootstrap: $(Child("run", CheckBox()))
""" end )

# ╔═╡ 7f720fc1-3004-41af-a570-46183ed4b646
md"""
The means and standard deviations of orbit parameters.
"""

# ╔═╡ 6926b12f-d5e9-475a-8248-2a4f0e0f18b5
md"""
If you want to see the distributions of parameters given by bootstrap similations, check the Appendix.
"""

# ╔═╡ 22ebdf59-a508-44f3-9db8-031b37c4446d
md"""
The cross-validation results are shown below.
"""

# ╔═╡ a40a805d-48a7-4451-b958-56ef041d3333
warning_box(md"This fitting method does not perform well for orbits with extremely long periods (>10000days). Please be cautious.")

# ╔═╡ 200c1f28-0eb3-459f-91f8-a0c38cb92bed
md"""
# Appendix
"""

# ╔═╡ 59069f4e-58fa-458a-8c47-4176e3970f43
md"""
The radial velocity with respect to phase plot, using the period with maximum power provided by lomb-scargle periodogram. The RV value for different instruments is randomly separated. 
"""

# ╔═╡ 5712b24d-cfff-4c95-8008-dfad4dc32a2c
md"""
RV data re-grouped by phasen is plotted below, using the period derived from the one planet keplerian model.
"""

# ╔═╡ 519b083c-aa8d-48da-b9ec-6e3cddf94d99
starid = searchsortedfirst(star_names,selected_star);

# ╔═╡ f3bbb76c-14cd-4566-935b-5c1f6949eaf5
begin
	star_name = star_names[starid]
	df_star = df_all |> @filter( _.Name == star_name ) |> DataFrame
end;

# ╔═╡ 2350b4a6-a538-430e-b832-4ecb4a458c4d
begin
	
	#df_star_by_inst=groupby(df_star,:Inst)
	#map(key->key.Inst, keys(df_star_by_inst))
	#combine(df_star_by_inst, nrow)
	df_star_by_inst = DataFrame()
	try
	df_star_by_inst = df_star |> @groupby( _.Inst ) |> @map( {bjd = _.d, rv = _.RVel, σrv = _.e_RVel, inst= key(_), nobs_inst=length(_) }) |> DataFrame;
	catch
	end
end;

# ╔═╡ 2de85bb1-881f-483a-bcdb-2109ed795ec5
 begin  # Make more useful observatory/instrument labels
	instrument_label = Dict(zip(["j","k","apf","lick"],["Keck (post)","Keck (pre)","APF","Lick"]))

	 #Instruments=[]
	#for (k, g) in pairs(df_star_by_inst)
    #     append!(Instruments,string(k.Inst))
    #end
	for k in keys(instrument_label)
		if k ∉ df_star_by_inst.inst
			delete!(instrument_label,k)
		end
	end
	instrument_label
end;

# ╔═╡ e3b46289-f602-4dee-81cc-91c8466dcd5a
df_star_by_inst |> @map( {:inst=>instrument_label[_.inst], :num_obs=>_.nobs_inst}) |> DataFrame

# ╔═╡ bfc1affc-b392-4884-a35f-55593af7db53
t_offset = 2455000; 

# ╔═╡ 6f000c6d-14a8-4d6d-97f8-d01f4dd516bc
begin
	plt_rv_all_inst = plot() #legend=:none, widen=true)
	local num_inst = size(df_star_by_inst,1)
	for inst in 1:num_inst
		rvoffset = mean(df_star_by_inst[inst,:rv])
		scatter!(plt_rv_all_inst,df_star_by_inst[inst,:bjd].-t_offset,
				df_star_by_inst[inst,:rv],
				yerr=collect(df_star_by_inst[inst,:σrv]),
				label=instrument_label[df_star_by_inst[inst,:inst]], markercolor=inst)
				#markersize=4*upscale, legendfontsize=upscale*12
	end
	xlabel!(plt_rv_all_inst,"Time (d)")
	ylabel!(plt_rv_all_inst,"RV (m/s)")
	title!(plt_rv_all_inst,"HD " * star_name )
	plt_rv_all_inst
end

# ╔═╡ 8f239f05-2610-42ce-92c6-968943418328
begin
	inst_list=[]

	for i in values(instrument_label)
		push!(inst_list,i)
	end
end

# ╔═╡ 383d9191-fc82-4f9c-81f5-837a67c71e9b
begin
	select_obs_cell_id = PlutoRunner.currently_running_cell_id[] |> string
	select_obs_cell_url = "#$(select_obs_cell_id)"
	md"""
Select first instrument's data to analyze: $(@bind inst1 Select(collect(values(instrument_label)); default=inst_list[1]))
"""
end

# ╔═╡ 10290d00-c667-46f3-a03e-8fb33842448a
if size(df_star_by_inst,1)>1 && ~fit_with_one 
	md"""
Select second instrument's data to analyze: $(@bind inst2 Select(collect(values(instrument_label)); default=inst_list[2]))
"""
end

# ╔═╡ 1ce61116-86e3-4488-b2ec-a0d7617bb4b3

md"""
Select the instrument for periodogram: $(@bind inst_per Select(collect(values(instrument_label)); default=inst_list[1])). Generate the Periodogram with a period up to $(@bind periodogram_range NumberField(0:0.1:500, default=50))% of the observation period of the selected instrument. 
"""

# ╔═╡ 32a6f1f8-c213-40c1-a6c1-5dc7ce87976e
begin
	inst_idx1 = findfirst(isequal(inst1),map(k->instrument_label[k], df_star_by_inst[:,:inst]));
	if @isdefined inst2 
		inst_idx2 = findfirst(isequal(inst2),map(k->instrument_label[k], df_star_by_inst[:,:inst]));
	end
	nothing
end

# ╔═╡ a435c976-de91-40a6-b536-9cf005027a3b
if size(df_star_by_inst,1)>0  # Warning: Assume exactly 2 instruments providing RV data
    data1 = (t=collect(df_star_by_inst[inst_idx1,:bjd]).-t_offset, rv=collect(df_star_by_inst[inst_idx1,:rv]), σrv=collect(df_star_by_inst[inst_idx1,:σrv]))
    if @isdefined inst_idx2
        data2 = (t=collect(df_star_by_inst[inst_idx2,:bjd]).-t_offset, rv=collect(df_star_by_inst[inst_idx2,:rv]), σrv=collect(df_star_by_inst[inst_idx2,:σrv]))
    else
        data2 = (t=Float64[], rv=Float64[], σrv=Float64[])
    end
    t_mean = (sum(data1.t)+sum(data2.t))/(length(data1.t).+length(data2.t))
    t_plt = range(minimum(vcat(data1.t,data2.t)), stop=maximum(vcat(data1.t,data2.t)), step=1.0)
end;

# ╔═╡ acf3b224-3aa7-4419-92a1-fd2b1ee582ce
md"""
For reference, the RMS of RV measurement uncertainty is $(round(sqrt(mean(vcat(data1.σrv,data2.σrv).^2)),sigdigits=3) ) m/s
"""

# ╔═╡ 963230d5-cef8-4a9c-aef4-319cbd75c207
if subtract_linear
begin
	plt_linear_subtracted = plot(xlabel="Time (d)", ylabel="RV (m/s)", title="RV data with linear results subtracted")
	scatter!(plt_linear_subtracted, data1.t,data1.rv,yerr=data1.σrv, label=instrument_label[df_star_by_inst[inst_idx1 ,:inst]], markercolor=inst_idx1)
	if @isdefined inst_idx2
		scatter!(plt_linear_subtracted, data2.t,data2.rv,yerr=data2.σrv, label=instrument_label[df_star_by_inst[inst_idx2 ,:inst]], markercolor=inst_idx2)
	end
	plt_linear_subtracted
end
end

# ╔═╡ d24fbe45-dc19-4946-9958-b3f26650d572
init_li=[0,mean(data1.rv),mean(data2.rv)];

# ╔═╡ 49a9a26e-6423-47c3-91ca-785c2ccafe24
inst_idx_per = findfirst(isequal(inst_per),map(k->instrument_label[k], df_star_by_inst[:,:inst]));

# ╔═╡ 1121402e-69d9-4d3b-bec9-c4011e226675
begin
	#s=df_star[(df_star[:,:Inst].=="k") .||(df_star[:,:Inst].=="j"),:RVel]
	#t=df_star[(df_star[:,:Inst].=="k").||(df_star[:,:Inst].=="j"),:d]
	#s=df_star[:,:RVel]
	#t=df_star[:,:d]
	#t=cat(data1.t,data2.t,dims=1)
	#s=cat(data1.rv,data2.rv,dims=1)
	#pgram_list=[]
	#maxPower_list=[]
	#for i =1:length(inst_list)
	#s=collect(df_star_by_inst[i,:rv])
	#t=collect(df_star_by_inst[i,:bjd])
	#D=maximum(t)-minimum(t) #Duration of observations
	#plan = LombScargle.plan(t, s,minimum_frequency=1/(D*periodogram_range/100))
	#pgram = lombscargle(plan)
	#push!(pgram_list,pgram)
	#push!(maxPower_list, LombScargle.findmaxpower(pgram))
	#end
	#id_pgram=1#findmax(maxPower_list)[2]
	#pgram=pgram_list[id_pgram]
	s=collect(df_star_by_inst[inst_idx_per,:rv])
	t=collect(df_star_by_inst[inst_idx_per,:bjd])
	D=maximum(t)-minimum(t) #Duration of observations
	plan = LombScargle.plan(t, s,minimum_frequency=1/(D*periodogram_range/100),normalization=:standard,oversampling=5)
	pgram = lombscargle(plan)
	plot(periodpower(pgram)...,label="Periodogram",xaxis=:log)
	fap=LombScargle.fapinv(pgram,0.0001)
	fap2=LombScargle.fapinv(pgram,0.01)
	hline!([fap],label="FAP=0.0001")
	hline!([fap2],label="FAP=0.01")
	title!("Periodogram")
	xlabel!("Period(d)")
	ylabel!("Power")
end

# ╔═╡ b51622f7-baf7-4ab5-861e-dd7f5a707fdb
if LombScargle.findmaxpower(pgram)<fap2
	warning_box(md"There is no period with FAP<0.01. That is, there is no obvious periodic signal. Be cautious in the following analysis.")
end

# ╔═╡ ec77269a-6cf8-420f-a7ef-75f15e30de28
md"""
Period with maximum power is $((LombScargle.findmaxperiod(pgram))). If you want to see how reasonable this period is, check the Appendix for the RV-phase diagram plotted with the period given by the periodogram. 
"""

# ╔═╡ 70e8d18a-fa92-4993-8994-c9a481141e1d
P0=LombScargle.findmaxperiod(pgram);

# ╔═╡ f46cbed8-dec6-4bb2-8942-4fe96ebea3e4
md"""
Set initial values. The initial guess for P is determined by the periodogram by default.\
P: $(@bind P_guess NumberField(1:0.0001:100000, default=P0[1]))



h: $(@bind h_guess NumberField(0:0.01:2π, default=0.01))
k: $(@bind k_guess NumberField(0:0.01:2π, default=0.01))
σj: $(@bind σj_guess NumberField(0:0.05:2π, default=3.0))
"""

# ╔═╡ 17b1760d-ed77-4773-8147-43245c70a962

begin
	
	plt_phase_all = plot(widen=false)
	num_inst = size(df_star_by_inst,1)
	for inst in 1:num_inst
		if length(df_star_by_inst[inst,:rv]) == 0 continue end
		rvoffset = mean(df_star_by_inst[inst,:rv]) .- 30 .* (inst-2)
		phase = mod.((df_star_by_inst[inst,:bjd].-t_offset)./P0,1.0)
		scatter!(plt_phase_all,phase,
				df_star_by_inst[inst,:rv].-rvoffset,
				yerr=collect(df_star_by_inst[inst,:σrv]),
				markersize=2, markerstrokewidth=0.5,
				label=instrument_label[df_star_by_inst[inst,:inst]])
	end
	#plot!(plt,t_plt,pred_1pl, label=:none)
	xlabel!(plt_phase_all,"Phase")
	ylabel!(plt_phase_all,"RV (m/s)")
	title!(plt_phase_all,"HD " * star_name *" with period from Periodogram" )
	plt_phase_all
end

# ╔═╡ 9e5a3e25-8f7a-469e-aa54-846beec8990f
# ╠═╡ disabled = true
#=╠═╡
begin
	res=vcat(resid1,resid2)
	s_res=res
	t_res=vcat(data1.t,data2.t)
	D_res=maximum(t_res)-minimum(t_res) #Duration of observations
	plan_res = LombScargle.plan(t_res, s_res,minimum_frequency=1/(D_res))
	pgram_res = lombscargle(plan_res)
end
  ╠═╡ =#

# ╔═╡ 769de9c4-6167-4550-be08-320c7c63fe3e
md"""
# Set Up
"""

# ╔═╡ 9d1eadf9-7fe9-4620-b866-4a442033a0b4
"""
`calc_mle_linear_model(A, y_obs, covar)`

Computes the maximum likelihood estimator for b for the linear model
`y = A b`
where measurements errors of `y_obs` are normally distributed and have covariance `covar`.
"""
function calc_mle_linear_model(A::AbstractMatrix, y_obs::AbstractVector, covar::AbstractMatrix)
	@assert size(A,1) == length(y_obs) == size(covar,1) == size(covar,2)
	@assert size(A,2) >= 1
	(A' * (covar \ A)) \ (A' * (covar \ y_obs) )
end;

# ╔═╡ 640c9851-73ba-47ed-84e4-6484f3b793b2
begin
	nobs1 = size(data1.t)[1]
	A1 = [ones(nobs1) data1.t-2455000*ones(nobs1)]
	b1=data1.rv
	covar1 = diagm(data1.σrv.^2)
	θ_mle1 = calc_mle_linear_model(A1, b1, covar1)

	nobs2 = size(data2.t)[1]
	A2 = [ones(nobs2) data2.t-2455000*ones(nobs2)]
	b2=data2.rv
	covar2 = diagm(data2.σrv.^2)
	θ_mle2 = calc_mle_linear_model(A2, b2, covar2)

	nothing
end

# ╔═╡ 59bc7b62-87f8-42bb-80d0-74d4cf2c9edf
begin
nobs=nobs1+nobs2;
nothing;
end

# ╔═╡ 5736df1b-a0c9-47b2-9cfe-faf01641ce26
function predict_linear_model(A::AbstractMatrix, b::AbstractVector)
	@assert size(A,2) == length(b)
	A*b
end;

# ╔═╡ 6cea0d3e-daed-4318-83b0-13bf9fe00d2a
function calc_true_anom(ecc_anom::Real, e::Real)
	true_anom = 2*atan(sqrt((1+e)/(1-e))*tan(ecc_anom/2))
end

# ╔═╡ c67d1588-36e4-44aa-a98b-8308bf57e1e0
"""
   `ecc_anom_init_guess_danby(mean_anomaly, eccentricity)`

Returns initial guess for the eccentric anomaly for use by iterative solvers of Kepler's equation for bound orbits.

Based on "The Solution of Kepler's Equations - Part Three"
Danby, J. M. A. (1987) Journal: Celestial Mechanics, Volume 40, Issue 3-4, pp. 303-312 (1987CeMec..40..303D)
"""
function ecc_anom_init_guess_danby(M::Real, ecc::Real)
	@assert -2π<= M <= 2π
	@assert 0 <= ecc <= 1.0
    if  M < zero(M)
		M += 2π
	end
    E = (M<π) ? M + 0.85*ecc : M - 0.85*ecc
end;

# ╔═╡ d892f267-8a1f-41e8-9104-573ec424abc9
"""
   `update_ecc_anom_laguerre(eccentric_anomaly_guess, mean_anomaly, eccentricity)`

Update the current guess for solution to Kepler's equation

Based on "An Improved Algorithm due to Laguerre for the Solution of Kepler's Equation"
   Conway, B. A.  (1986) Celestial Mechanics, Volume 39, Issue 2, pp.199-211 (1986CeMec..39..199C)
"""
function update_ecc_anom_laguerre(E::Real, M::Real, ecc::Real)
  #es = ecc*sin(E)
  #ec = ecc*cos(E)
  (es, ec) = ecc .* sincos(E)  # Does combining them provide any speed benefit?
  F = (E-es)-M
  Fp = one(M)-ec
  Fpp = es
  n = 5
  root = sqrt(abs((n-1)*((n-1)*Fp*Fp-n*F*Fpp)))
  denom = Fp>zero(E) ? Fp+root : Fp-root
  return E-n*F/denom
end;

# ╔═╡ 5dcd30d7-dd60-4b15-96a6-4b222e55d779
begin
	calc_ecc_anom_cell_id = PlutoRunner.currently_running_cell_id[] |> string
	calc_ecc_anom_url = "#$(calc_ecc_anom_cell_id)"
	"""
	   `calc_ecc_anom( mean_anomaly, eccentricity )`
	   `calc_ecc_anom( param::Vector )`

	Estimates eccentric anomaly for given 'mean_anomaly' and 'eccentricity'.
	If passed a parameter vector, param[1] = mean_anomaly and param[2] = eccentricity.

	Optional parameter `tol` specifies tolerance (default 1e-8)
	"""
	function calc_ecc_anom end

	function calc_ecc_anom(mean_anom::Real, ecc::Real; tol::Real = 1.0e-8)
	  	if !(0 <= ecc <= 1.0)
			println("mean_anom = ",mean_anom,"  ecc = ",ecc)
		end
		@assert 0 <= ecc <= 1.0
		@assert 1e-16 <= tol < 1
	  	M = rem2pi(mean_anom,RoundNearest)
	    E = ecc_anom_init_guess_danby(M,ecc)
		local E_old
	    max_its_laguerre = 200
	    for i in 1:max_its_laguerre
	       E_old = E
	       E = update_ecc_anom_laguerre(E_old, M, ecc)
	       if abs(E-E_old) < tol break end
	    end
	    return E
	end

	function calc_ecc_anom(param::Vector; tol::Real = 1.0e-8)
		@assert length(param) == 2
		calc_ecc_anom(param[1], param[2], tol=tol)
	end;
end

# ╔═╡ 0f2ee09b-1d31-4a34-9b0a-22e8b91f4303
begin
	""" Calculate RV from t, P, K, e, ω and M0	"""
	function calc_rv_keplerian end
	calc_rv_keplerian(t, p::Vector) = calc_rv_keplerian(t, p...)
	function calc_rv_keplerian(t, P,K,e,ω,M0)
		mean_anom = t*2π/P-M0
		ecc_anom = calc_ecc_anom(mean_anom,e)
		true_anom = calc_true_anom(ecc_anom,e)
		rv = K/sqrt((1-e)*(1+e))*(cos(ω+true_anom)+e*cos(ω))
	end
end

# ╔═╡ 74170a94-2014-45fd-9f1d-d476b1febf14
begin
	""" Calculate RV from t, P, K, e, ω, M0	and C"""
	function calc_rv_keplerian_plus_const end
	calc_rv_keplerian_plus_const(t, p::Vector) = calc_rv_keplerian_plus_const(t, p...)

	function calc_rv_keplerian_plus_const(t, P,K,e,ω,M0,C)
		calc_rv_keplerian(t, P,K,e,ω,M0) + C
	end
end

# ╔═╡ 2293d249-407a-4578-8eb8-377a3ef44e68
function loss_li(θ)
		(slope, C1,C2) = θ
		
		C2=C1
		rv_model1 = slope.*data1.t.+C1
		loss =sum( (rv_model1.-data1.rv).^2 ./ data1.σrv.^2 )
	if @isdefined inst_idx2	
	rv_model2 = data2.t*slope.+C2
		loss += sum( (rv_model2.-data2.rv).^2 ./ data2.σrv.^2  )
	end
		return loss
	end

# ╔═╡ 75e5df52-18ab-4046-a2a0-fdaea68d9543
result_li=Optim.optimize(loss_li, init_li, BFGS(), autodiff=:forward);

# ╔═╡ 7e9f805c-80c4-472f-a3ea-f73a732b8a57
md"""
The loss function value with best-fit parameters is $(loss_li(result_li.minimizer)). The number of observations is $(nobs1+nobs2).
"""

# ╔═╡ 36e7d4ab-735e-4517-bf2e-ae1db69d227e
loss_linear = loss_li(result_li.minimizer);

# ╔═╡ d63c6ac8-3623-426e-b25f-668ffce47ddd
let
	expected_distribution = Chisq(nobs1+nobs2)
	est_max_y = maximum(pdf.(expected_distribution,range(nobs1+nobs2/2,stop=1.5*(nobs1+nobs2),length=100)))

	plt = plot(expected_distribution, xlabel="loss", ylabel="Probability", legend=:none)
	plot!(plt,[loss_linear,loss_linear],[0,est_max_y], linestyle=:dot, linecolor=:black)
end

# ╔═╡ 44935b70-11a6-4804-a7cb-150d8660441b
if loss_linear > nobs*5
	warning_box(md"Linear model is probably not a good fit to the selected datasets.")
end

# ╔═╡ 3e653573-04db-4913-8771-2c23fe0e01a1
begin
	slope_li=result_li.minimizer[1];
	C1_li=result_li.minimizer[2];
	C2_li=result_li.minimizer[3];
	nothing
end

# ╔═╡ b5fc8991-ae2d-4db4-a7f9-5da11a21e094
begin
	plt_linear = plot(xlabel="Time (d)", ylabel="RV (m/s)", title="RV data with linear model")
	scatter!(plt_linear, data1.t,data1.rv.-C1_li,yerr=data1.σrv, label=instrument_label[df_star_by_inst[inst_idx1 ,:inst]], markercolor=inst_idx1)
	if @isdefined inst_idx2
		scatter!(plt_linear, data2.t,data2.rv.-C1_li,yerr=data2.σrv, label=instrument_label[df_star_by_inst[inst_idx2 ,:inst]], markercolor=inst_idx2)
	
	plot!([minimum([minimum(data2.t),minimum(data1.t)]),maximum([maximum(data2.t),maximum(data1.t)])],[slope_li*minimum([minimum(data2.t),minimum(data1.t)]),slope_li*maximum([maximum(data2.t),maximum(data1.t)])],label="Linear fit")
	plt_linear
	else
	plot!([minimum(data1.t),maximum(data1.t)],[slope_li*minimum(data1.t),slope_li*maximum(data1.t)],label="Linear fit")
	plt_linear
	end
end

# ╔═╡ 7f8751c3-b087-401e-a098-25460bf496d9
begin
	resid_li1=data1.rv.-C1_li-slope_li.*data1.t
	resid_li2=data2.rv.-C1_li-slope_li.*data2.t
	
	plt_resid_li = plot(legend=:none, widen=true)
	scatter!(plt_resid_li,data1.t,
				resid_li1,
				yerr=data1.σrv, markercolor=inst_idx1)
	
	xlabel!(plt_resid_li,"Time (d)")
	ylabel!(plt_resid_li,"RV (m/s)")

	

	if @isdefined inst_idx2
	scatter!(plt_resid_li,data2.t,
				resid_li2,
				yerr=data2.σrv, markercolor=inst_idx2)
	
	end
	
	title!(plt_resid_li,"HD " * star_name * " (residuals to linear model)")
	plt_resid_li
end

# ╔═╡ c2c32163-dcbc-4bfc-a172-e3d750ce6c4b
if subtract_linear
	for i=1:length(data1.rv)
		data1.rv[i]=data1.rv[i]-C1_li-slope_li*data1.t[i]
	end
	for i=1:length(data2.rv)
		data2.rv[i]=data2.rv[i]-C2_li-slope_li*data2.t[i]
	end
end

# ╔═╡ 0372f0ee-b238-403e-b045-1ade617a271d
""" Calculate RV from t, P, K, e, ω, M0	and C with optional slope and t_mean"""
function model_1pl(t, P, K, e, ω, M, C; slope=0.0, t_mean = 0.0)
    calc_rv_keplerian(t-t_mean,P,K,e,ω,M) + C + slope * (t-t_mean)
end


# ╔═╡ 61379c75-6969-4652-b6ce-7e1e992f1f23
""" Convert vector of (P,K,h,k,ω+M0) to vector of (P, K, e, ω, M0) """
function PKhkωpM_to_PKeωM(x::Vector) 
    (P, K, h, k, ωpM) = x
    ω = atan(h,k)
    return [P, K, sqrt(h^2+k^2), ω, ωpM-ω]
end

# ╔═╡ 1d397e31-740e-41e7-ab4c-d6009752ee33
""" Convert vector of (P,K,h,k,M0-ω) to vector of (P, K, e, ω, M0) """
function PKhkωMmω_to_PKeωM(x::Vector)
	(P, K, h, k, ωmM) = x
	e = sqrt(h^2+k^2)
	ω = atan(h,k)
	return [P, K, e, ω, ωmM+ω]
end

# ╔═╡ 71e66a2a-be4d-48f3-a868-bcd8f7647e22
function make_loss_1pl(data1,data2; t_mean=0)
	function loss_1pl(θ)
		(P1, K1, h1, k1, Mpω1, C1,C2, σj ) = θ
		( P1, K1, e1, ω1, M1 ) = PKhkωMmω_to_PKeωM([P1, K1, h1, k1, Mpω1])
		if e1>1 return 1e6*e1 end
		rv_model1 = model_1pl.(data1.t,P1,K1,e1,ω1,M1,C1,t_mean=t_mean)
		loss = 0.5*sum( (rv_model1.-data1.rv).^2 ./ (data1.σrv.^2 .+ σj^2) )
		rv_model2 = model_1pl.(data2.t,P1,K1,e1,ω1,M1,C2, t_mean=t_mean)
		loss += 0.5*sum( (rv_model2.-data2.rv).^2 ./ (data2.σrv.^2 .+ σj^2) )
		loss += 0.5*sum(log.(2π*(data1.σrv.^2 .+σj^2)))
		loss += 0.5*sum(log.(2π*(data2.σrv.^2 .+σj^2)))
		return loss
	end
	#function loss_1pl(θ) 
    #(P1, K1, h1, k1, Mpω1, C1, C2, slope, σj ) = θ
    #( P1, K1, e1, ω1, M1 ) = PKhkωpM_to_PKeωM([P1, K1, h1, k1, Mpω1])
    #if e1>1 return 1e6*e1 end
    #rv_model1 = model_1pl.(data1.t,P1,K1,e1,ω1,M1,C1, slope=slope, t_mean=t_mean)
    #loss = 0.5*sum(((rv_model1.-data1.rv)./(data1.σrv.+σj^2)).^2)
    #rv_model2 = model_1pl.(data2.t,P1,K1,e1,ω1,M1,C2, slope=slope, t_mean=t_mean)
    #loss += 0.5*sum((rv_model2.-data2.rv).^2 ./(data2.σrv.^2 .+σj^2))
    #loss += 0.5*sum(log.(2π*(data1.σrv.^2 .+σj^2)))
    #loss += 0.5*sum(log.(2π*(data2.σrv.^2 .+σj^2)))
    #return loss
#end

end

# ╔═╡ 514730e4-f589-4841-a4db-4b39e746e6a9
function find_best_1pl_fit(θinit::AbstractVector, loss::Function; num_init_phases::Integer=1, num_init_ωs::Integer=1, f_abstol::Real = 1e-2 )
	@assert 1 <= num_init_phases <= 32
	@assert 1 <= num_init_ωs <= 8
	result_list = Array{Any}(undef,num_init_phases, num_init_ωs)
	θinit_list = fill(θinit,num_init_phases, num_init_ωs)
	e_base = sqrt(θinit[3]^2+θinit[4]^2)
	ω_base = atan(θinit[3],θinit[4])
	for i in 1:num_init_phases
		for j in 1:num_init_ωs
		Δω = (j-1)/num_init_ωs * 2π
		θinit_list[i,j][3] = e_base*sin(ω_base + Δω)
		θinit_list[i,j][4] = e_base*cos(ω_base + Δω)
		θinit_list[i,j][5] += (i-1)/num_init_phases * 2π - Δω
		θinit_list[i,j][5] = mod(θinit_list[i,j][5],2π)
		try
			result_list[i,j] = Optim.optimize(loss ,θinit_list[i,j], BFGS(), autodiff=:forward, Optim.Options(f_abstol=f_abstol));
		catch
			result_list[i,j] = (;minimum=Inf)
		end
		end
	end
	best_result_id = argmin(map(r->r.minimum, result_list))
	result = result_list[best_result_id]
end

# ╔═╡ 96851df0-0aad-44ac-830a-b5062917f4cf
function fit_general_linear_least_squares( design_mat::ADM, covar_mat::APD, obs::AA ) where { ADM<:AbstractMatrix, APD<:AbstractPDMat, AA<:AbstractArray }
   Xt_inv_covar_X = Xt_invA_X(covar_mat,design_mat)
   X_inv_covar_y =  design_mat' * (covar_mat \ obs)
   AB_hat = Xt_inv_covar_X \ X_inv_covar_y                   # standard GLS
end

# ╔═╡ 167cc8ad-cb04-47a1-bb25-08e9c345b24e
function calc_design_matrix_circ!(result::AM, period, times::AV) where { R1<:Real, AM<:AbstractMatrix{R1}, AV<:AbstractVector{R1} }
        n = length(times)
        @assert size(result) == (n, 2)
        for i in 1:n
                ( result[i,1], result[i,2] ) = sincos(2π/period .* times[i])
        end
        return result
end

# ╔═╡ 5edca12e-e708-43d3-b0dc-4b70c1c4ea70
function calc_design_matrix_circ(period, times::AV) where { R1<:Real, AV<:AbstractVector{R1} }
        n = length(times)
        dm = zeros(n,2)
        calc_design_matrix_circ!(dm,period,times)
        return dm
end

# ╔═╡ acff1e9d-038c-4267-9f85-37e21122988c
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	param_fit_linear = fit_general_linear_least_squares(
       calc_design_matrix_circ(P_guess,data1.t), PDiagMat(data1.σrv), data1.rv)
	K_guess_def = sqrt(param_fit_linear[1]^2+param_fit_linear[2]^2)
	phase_guess = atan(param_fit_linear[1],param_fit_linear[2]);
	nothing
end
  ╠═╡ =#

# ╔═╡ aae7b265-1ff7-4d5b-8d93-f0bb0bce575a
#=╠═╡
md"""
Sometimes the estimation of K can be very sensitive to P. The default K is given by fitting data from the first chosen dataset to a circular model. If the default K looks suspicious or gives bad results, you can change K as you wish. 
K: $(@bind K_guess NumberField(0:0.01:200, default=K_guess_def))
"""
  ╠═╡ =#

# ╔═╡ 58b75f36-b9f5-4192-a796-73f65e497df2
#=╠═╡
begin
	C1_guess=mean(data1.rv)
	C2_guess=mean(data2.rv)
	#init_guess=[4.23079,54.0364,-0.01,0.01,4.79533,C1_guess,C2_guess,3.0]
	init_guess = [P_guess, K_guess, h_guess, k_guess, mod(phase_guess-atan(h_guess, k_guess),2π), C1_guess,C2_guess, σj_guess]
end
  ╠═╡ =#

# ╔═╡ 34fa0f41-fae2-4854-8055-d9ee476c3eef
#=╠═╡
	result = find_best_1pl_fit(init_guess, make_loss_1pl(data1,data2, t_mean=t_mean), num_init_phases = 1, num_init_ωs=4);
  ╠═╡ =#

# ╔═╡ db3dca2a-86ed-4c6e-9477-175092eb538f
#=╠═╡
result.minimizer
  ╠═╡ =#

# ╔═╡ 255ce271-2011-4082-acaa-bfd4d972dc7b
#=╠═╡
if @isdefined result
let
	
	#upscale
	plt = plot() #legend=:none, widen=true)
	num_inst = size(df_star_by_inst,1)
	rvoffset = zeros(4) # [result2.minimizer[11], result2.minimizer[12], 0, 0]
	rvoffset[inst_idx1]= result.minimizer[6]
	
	
		
		#plot!(t_plt,model_1pl)
		#plot!(t_plt.-2450000,model_2pl)
		#scatter!(df_star_by_inst[2,:bjd].-2450000,df_star_by_inst[2,:rv], yerr=collect(df_star_by_inst[2,:σrv]),label=:none)
		#scatter!(df_star_by_inst[3,:bjd].-2450000,df_star_by_inst[3,:rv], yerr=collect(df_star_by_inst[3,:σrv]),label=:none)
		#scatter!(df_star_by_inst[4,:bjd].-2450000,df_star_by_inst[4,:rv], yerr=collect(df_star_by_inst[4,:σrv]),label=:none)
	
	pred_1pl = map(t->model_1pl(t,PKhkωMmω_to_PKeωM(result.minimizer[1:5])...,0.0,slope=0.0, t_mean=t_mean),t_plt);
	plot!(t_plt,pred_1pl, label=:none)
	xlabel!("Time (d)")
	ylabel!("RV (m/s)")
	title!("HD " * star_name )
	#savefig(plt,joinpath(homedir(),"Downloads","RvEx.pdf"))
	if @isdefined inst_idx2
		rvoffset[inst_idx2]= result.minimizer[7]
		scatter!(plt, data2.t,data2.rv.-rvoffset[inst_idx2],yerr=data2.σrv, label=instrument_label[df_star_by_inst[inst_idx2 ,:inst]], markercolor=inst_idx2)
	end
	slope = 0 #result.minimizer[8]
	scatter!(plt, data1.t,data1.rv.-rvoffset[inst_idx1],yerr=data1.σrv, label=instrument_label[df_star_by_inst[inst_idx1 ,:inst]], markercolor=inst_idx1)
	plt
end
	
end
  ╠═╡ =#

# ╔═╡ 34304691-cdd8-4cc0-a33e-e166c434eb4d
#=╠═╡
begin
	plt_phase_all_after = plot(widen=false)
	for inst in 1:num_inst
		if length(df_star_by_inst[inst,:rv]) == 0 continue end
		rvoffset = mean(df_star_by_inst[inst,:rv]) .- 30 .* (inst-2)
		phase = mod.((df_star_by_inst[inst,:bjd].-t_offset)./result.minimizer[1],1.0)
		scatter!(plt_phase_all_after,phase,
				df_star_by_inst[inst,:rv].-rvoffset,
				yerr=collect(df_star_by_inst[inst,:σrv]),
				markersize=2, markerstrokewidth=0.5,
				label=instrument_label[df_star_by_inst[inst,:inst]])
	end
	#plot!(plt,t_plt,pred_1pl, label=:none)
	xlabel!(plt_phase_all_after,"Phase")
	ylabel!(plt_phase_all_after,"RV (m/s)")
	title!(plt_phase_all_after,"HD " * star_name * " with the fitted period." )
	plt_phase_all_after
end
  ╠═╡ =#

# ╔═╡ 4a3292ab-cfa2-4b58-ba66-19a23f8f2a23
#=╠═╡
begin
pred_1pl1 = map(t->model_1pl(t,PKhkωMmω_to_PKeωM(result.minimizer[1:5])...,0.0, t_mean=t_mean),data1.t)
			pred_1pl2 = map(t->model_1pl(t,PKhkωMmω_to_PKeωM(result.minimizer[1:5])...,0.0, t_mean=t_mean),data2.t)
			resid1 = data1.rv.-result.minimizer[6]-pred_1pl1
			resid2 = data2.rv.-result.minimizer[7]-pred_1pl2
nothing
end
  ╠═╡ =#

# ╔═╡ 3a51761b-8659-4c42-9612-28d6daaf7404
#=╠═╡
if @isdefined resid1
	plt_resid = plot(legend=:none, widen=true)
	scatter!(plt_resid,data1.t,
				resid1,
				yerr=data1.σrv, markercolor=inst_idx1)
	
	xlabel!(plt_resid,"Time (d)")
	ylabel!(plt_resid,"RV (m/s)")

	#plt_resid_phase = plot(legend=:none, widen=false)
	#phase = mod.(data1.t ./ result.minimizer[1],1.0)
	#scatter!(plt_resid_phase,phase,
				#resid1,
				#yerr=data1.σrv, markercolor=inst_idx1)

	if @isdefined inst_idx2
	scatter!(plt_resid,data2.t,
				resid2,
				yerr=data2.σrv, markercolor=inst_idx2)
	#scatter!(plt_resid_phase,phase,
				#resid2,
				#yerr=data2.σrv, markercolor=inst_idx2)
	end
	#xlabel!(plt_resid_phase,"Phase")
	#ylabel!(plt_resid_phase,"RV (m/s)")
	title!(plt_resid,"HD " * star_name * " (residuals to 1 planet model)")
	#plot(plt_resid, plt_resid_phase, layout=(2,1) )
	plot(plt_resid)
end
  ╠═╡ =#

# ╔═╡ e8222797-77df-40d2-b9bd-279ed3c12cde
#=╠═╡
md"""
The mean measurement uncertainty is $(round(mean(vcat(data1.σrv,data2.σrv)),sigdigits=3)) m/s. The RMS RV residual is $(round(sqrt(mean(vcat(resid1,resid2).^2)),sigdigits=3)) m/s.
"""
  ╠═╡ =#

# ╔═╡ 06eede8c-1dda-4a15-a69c-6019ee11269d
#=╠═╡
if sqrt(mean(vcat(resid1,resid2).^2))<2*mean(vcat(data1.σrv,data2.σrv)) 
	correct("The model seems to be a good fit!")
elseif  8*mean(vcat(data1.σrv,data2.σrv)) >sqrt(mean(vcat(resid1,resid2).^2)) > 2*mean(vcat(data1.σrv,data2.σrv)) 
	warning_box(md"1 planet model is probably not a good fit to the selected datasets. Please check the fitting results. Please also check the residuals to see if there is another planet.")
elseif 8*mean(vcat(data1.σrv,data2.σrv))<sqrt(mean(vcat(resid1,resid2).^2))
	warning_box(md"1 planet model is not a good fit to the selected datasets. Please check the fitting results. Please also check the residuals to see if there is another planet.")
end
  ╠═╡ =#

# ╔═╡ 7d7fba67-c78a-44fe-bd9b-d8d4967b32c7
#=╠═╡
if bootstrap.run
	num_bootstrap_samples = bootstrap.num_bootstrap_samples
results_bootstrap = Array{Any}(undef,num_bootstrap_samples)
	rms_bootstrap = zeros(num_bootstrap_samples)
	rms_bootstrap_train = zeros(num_bootstrap_samples)
	testset_length = zeros(num_bootstrap_samples)
	global not_converge_flag=0
	not_converge_index=[]
	for i = 1:num_bootstrap_samples
		# Select sample of points to fit
		idx1 = sample(1:length(data1.t),length(data1.t))
		idx2 = sample(1:length(data2.t),length(data2.t))
		# Create NamedTuple with view of resampled data
		data_tmp1 = (;t=view(data1.t,idx1), rv=view(data1.rv,idx1), σrv=view(data1.σrv,idx1))
		data_tmp2 = (;t=view(data2.t,idx2), rv=view(data2.rv,idx2), σrv=view(data2.σrv,idx2))
	
		loss_tmp = make_loss_1pl(data_tmp1,data_tmp2, t_mean=t_mean)
		# Attempt to find best-fit parameters results for resampled data
		try 
			find_best_1pl_fit(result.minimizer[1:8], loss_tmp, num_init_phases=1, num_init_ωs=4).minimizer
		catch e
			if e==ErrorException("type NamedTuple has no field minimizer")
				not_converge_flag=1
				println("Does not converge")
				push!(not_converge_index,i)
				continue
			end

		end
		results_bootstrap[i]=find_best_1pl_fit(result.minimizer[1:8], loss_tmp, num_init_phases=1, num_init_ωs=4)
		# Evaluate residuals for cross validation
		if hasfield(typeof(results_bootstrap[i]),:minimizer)
			idx_train1 = filter(i->(i∈idx1), 1:length(data1.t))
			idx_train2 = filter(i->(i∈idx2), 1:length(data2.t))
			if length(idx_train1) == 0 continue end
			local pred_1pl1 = map(t->model_1pl(t,PKhkωMmω_to_PKeωM(results_bootstrap[i].minimizer[1:5])...,0.0, t_mean=t_mean),view(data1.t,idx_train1))
			local pred_1pl2 = map(t->model_1pl(t,PKhkωMmω_to_PKeωM(results_bootstrap[i].minimizer[1:5])...,0.0, t_mean=t_mean),view(data2.t,idx_train2))
			local resid1 = view(data1.rv,idx_train1).-results_bootstrap[i].minimizer[6]-pred_1pl1
			local resid2 = view(data2.rv,idx_train2).-results_bootstrap[i].minimizer[7]-pred_1pl2
			resid=vcat(resid1,resid2)
			rms_bootstrap_train[i] = sqrt(mean(resid.^2))
				
			idx_test1 = filter(i->!(i∈idx1), 1:length(data1.t))
			idx_test2 = filter(i->!(i∈idx2), 1:length(data2.t))
			testset_length[i] = length(idx_test1)+length(idx_test2)
			if testset_length[i] == 0 continue end
			local pred_1pl1 = map(t->model_1pl(t,PKhkωMmω_to_PKeωM(results_bootstrap[i].minimizer[1:5])...,0.0, t_mean=t_mean),view(data1.t,idx_test1))
			local pred_1pl2 = map(t->model_1pl(t,PKhkωMmω_to_PKeωM(results_bootstrap[i].minimizer[1:5])...,0.0, t_mean=t_mean),view(data2.t,idx_test2))
			local resid1 = view(data1.rv,idx_test1).-results_bootstrap[i].minimizer[6]-pred_1pl1
			local resid2 = view(data2.rv,idx_test2).-results_bootstrap[i].minimizer[7]-pred_1pl2
			resid=vcat(resid1,resid2)
			rms_bootstrap[i] = sqrt(mean(resid.^2))
			
		
	end
	end
	results_bootstrap=results_bootstrap[Not(not_converge_index)]
	nothing
end
  ╠═╡ =#

# ╔═╡ 2d40a7b2-f650-45e0-978e-eca24954faa6
#=╠═╡
if @isdefined rms_bootstrap
	plt_bootstrap_resid = plot(xlabel="RMS RV Residuals (m/s)", ylabel="Samples")
	histogram!(plt_bootstrap_resid, rms_bootstrap, nbins=40, label="Test points", alpha=0.5)
	histogram!(plt_bootstrap_resid, rms_bootstrap_train, nbins=40, label="Training points", alpha=0.5)
end
  ╠═╡ =#

# ╔═╡ 79e5a375-3c16-4443-86b9-d3bbad6101d4
#=╠═╡
if @isdefined results_bootstrap
	plt_title = plot(title = "Bootstrap Results", grid = false, showaxis = false, ticks=:none, bottom_margin = -25Plots.px)

	local Psample = map(r->r.minimizer[1],results_bootstrap)
	P_mean_bootstrap = mean(Psample)
	P_std_bootstrap = std(Psample)
	plt_P_hist = plot(xlabel="P (d)",ylabel="Samples",xticks=
	optimize_ticks(minimum(Psample),maximum(Psample),k_max=3)[1])
	histogram!(plt_P_hist,Psample, label=:none, nbins=50)

	local Ksample = map(r->r.minimizer[2],results_bootstrap)
	K_mean_bootstrap = mean(Ksample)
	K_std_bootstrap = std(Ksample)
	plt_K_hist = plot(xlabel="K (m/s)", ylabel="Samples",xticks=
	optimize_ticks(minimum(Ksample),maximum(Ksample),k_max=3)[1])
	histogram!(plt_K_hist,Ksample, label=:none, nbins=50)

	local esample = map(r->PKhkωMmω_to_PKeωM(r.minimizer)[3],results_bootstrap)
	e_mean_bootstrap = mean(esample)
	e_std_bootstrap = std(esample)
	plt_e_hist = plot(xlabel="e", ylabel="Samples",xticks=
	optimize_ticks(minimum(esample),maximum(esample),k_max=3)[1])
	histogram!(plt_e_hist,esample, label=:none, nbins=50)

	local ωsample = map(r->PKhkωMmω_to_PKeωM(r.minimizer)[4],results_bootstrap)
	ω_mean_bootstrap = mean(ωsample)
	ω_std_bootstrap = std(ωsample)
	plt_ω_hist = plot(xlabel="ω", ylabel="Samples",xticks=
	optimize_ticks(minimum(ωsample),maximum(ωsample),k_max=3)[1])
	histogram!(plt_ω_hist,ωsample, label=:none, nbins=50)

	h_mean_bootstrap = mean(esample.*sin.(ωsample))
	h_std_bootstrap = std(esample.*sin.(ωsample))
	k_mean_bootstrap = mean(esample.*cos.(ωsample))
	k_std_bootstrap = std(esample.*cos.(ωsample))

	local Mmωsample = map(r->PKhkωMmω_to_PKeωM(r.minimizer)[5],results_bootstrap)
	Mmω_mean_bootstrap = mean(Mmωsample)
	Mmω_std_bootstrap = std(Mmωsample)
	plt_Mmω_hist = plot(xlabel="M₀-ω", ylabel="Samples",xticks=
	optimize_ticks(minimum(Mmωsample),maximum(Mmωsample),k_max=3)[1])
	histogram!(plt_Mmω_hist,ωsample, label=:none, nbins=50)

	local C1sample = map(r->r.minimizer[6],results_bootstrap)
	C1_mean_bootstrap = mean(C1sample)
	C1_std_bootstrap = std(C1sample)
	plt_C1_hist = plot(xlabel="C1 (m/s)", ylabel="Samples",xticks=
	optimize_ticks(minimum(C1sample),maximum(C1sample),k_max=2)[1])
	histogram!(plt_C1_hist,C1sample, label=:none, nbins=50)

	if @isdefined inst_idx2
		local C2sample = map(r->r.minimizer[7],results_bootstrap)
		C2_mean_bootstrap = mean(C2sample)
		C2_std_bootstrap = std(C2sample)
		plt_C2_hist = plot(xlabel="C2 (m/s)", ylabel="Samples",xticks=
		optimize_ticks(minimum(C2sample),maximum(C2sample),k_max=2)[1])
		histogram!(plt_C2_hist,C2sample, label=:none, nbins=50)
	end

	local σjsample = map(r->r.minimizer[8],results_bootstrap)
	σj_mean_bootstrap = mean(σjsample)
	σj_std_bootstrap = std(σjsample)
	plt_σj_hist = plot(xlabel="σⱼ (m/s)", ylabel="Samples",xticks=
	optimize_ticks(minimum(σjsample),maximum(σjsample),k_max=3)[1])
	histogram!(plt_σj_hist,σjsample, label=:none, nbins=50)
	if @isdefined inst_idx2
		bootstrap_results_df = DataFrame(:parameter=>[:P, :K, :h, :k, :e, :ω, :Mmω, :C1,:C2, :σj], :mean=>[P_mean_bootstrap, K_mean_bootstrap, h_mean_bootstrap, k_mean_bootstrap, e_mean_bootstrap, ω_mean_bootstrap, Mmω_mean_bootstrap, C1_mean_bootstrap, C2_mean_bootstrap, σj_mean_bootstrap], :std=>[P_std_bootstrap,K_std_bootstrap,h_std_bootstrap,k_std_bootstrap,e_std_bootstrap,ω_std_bootstrap,Mmω_std_bootstrap,C1_std_bootstrap, C2_std_bootstrap,σj_std_bootstrap])
	
		plot(plt_title, plt_P_hist, plt_K_hist, plt_e_hist, plt_ω_hist, plt_C1_hist, plt_C2_hist, plt_σj_hist, layout=@layout([A{0.01h}; [B C; D E; F G;H ]]), size=(600,600) )
	else
		bootstrap_results_df = DataFrame(:parameter=>[:P, :K, :h, :k, :e, :ω, :Mmω, :C1, :σj], :mean=>[P_mean_bootstrap, K_mean_bootstrap, h_mean_bootstrap, k_mean_bootstrap, e_mean_bootstrap, ω_mean_bootstrap, Mmω_mean_bootstrap, C1_mean_bootstrap,  σj_mean_bootstrap], :std=>[P_std_bootstrap,K_std_bootstrap,h_std_bootstrap,k_std_bootstrap,e_std_bootstrap,ω_std_bootstrap,Mmω_std_bootstrap,C1_std_bootstrap, σj_std_bootstrap])

		plot(plt_title, plt_P_hist, plt_K_hist, plt_e_hist, plt_ω_hist, plt_C1_hist,  plt_σj_hist, layout=@layout([A{0.01h}; [B C; D E; F G ]]), size=(600,600) )
	end
end
  ╠═╡ =#

# ╔═╡ c344a8d1-da07-4206-87d2-92924e133001
#=╠═╡
if @isdefined bootstrap_results_df
   bootstrap_results_df
end
  ╠═╡ =#

# ╔═╡ e30755c3-e701-4fd6-a244-6755ed6603d3
#=╠═╡
if @isdefined not_converge_flag
	if not_converge_flag==1
		warning_box(md"Several samples do not converge.")
	end
end
  ╠═╡ =#

# ╔═╡ 66446e92-0d84-4d2b-9804-a563fbd8e30b
# ╠═╡ disabled = true
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	plt_1inst = plot(xlabel="Time (d)", ylabel="RV (m/s)", title="RVs from instrument being fit")
	scatter!(plt_1inst, df_star_by_inst[inst_idx1,:bjd].-t_offset,df_star_by_inst[inst_idx1,:rv],yerr=df_star_by_inst[inst_idx1,:σrv], label=instrument_label[df_star_by_inst[inst_idx1 ,:inst]], markercolor=inst_idx1)
	if @isdefined inst_idx2
		scatter!(plt_1inst, df_star_by_inst[inst_idx2,:bjd].-t_offset,df_star_by_inst[inst_idx2,:rv],yerr=df_star_by_inst[inst_idx2,:σrv], label=instrument_label[df_star_by_inst[inst_idx2 ,:inst]], markercolor=inst_idx2)
	end
end
  ╠═╡ =#

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
Downloads = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
LombScargle = "fc60dff9-86e7-5f2f-a8a0-edeadbb75bd9"
MCMCChains = "c7f686f2-ff18-58e9-bc7b-31028e88f75d"
MarkdownLiteral = "736d6165-7244-6769-4267-6b50796e6954"
Optim = "429524aa-4258-5aef-a3af-852621145aeb"
PDMats = "90014a1f-27ba-587c-ab20-58faa44d9150"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Query = "1a8c2f83-1ff3-5112-b086-8aa67b057ba1"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
StatsPlots = "f3b207a7-027a-5e70-b257-86293d7955fd"
Turing = "fce5fe82-541a-59a6-adf8-730c64b5f9a0"

[compat]
CSV = "~0.10.7"
DataFrames = "~1.4.4"
Distributions = "~0.25.79"
ForwardDiff = "~0.10.32"
LaTeXStrings = "~1.3.0"
LombScargle = "~1.0.3"
MCMCChains = "~5.5.0"
MarkdownLiteral = "~0.1.1"
Optim = "~1.7.4"
PDMats = "~0.11.16"
Plots = "~1.37.1"
PlutoTeachingTools = "~0.2.5"
PlutoUI = "~0.7.49"
Query = "~1.0.0"
StatsBase = "~0.33.21"
StatsPlots = "~0.15.4"
Turing = "~0.22.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.3"
manifest_format = "2.0"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "69f7020bd72f069c219b5e8c236c1fa90d2cb409"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.2.1"

[[deps.AbstractMCMC]]
deps = ["BangBang", "ConsoleProgressMonitor", "Distributed", "Logging", "LoggingExtras", "ProgressLogging", "Random", "StatsBase", "TerminalLoggers", "Transducers"]
git-tree-sha1 = "5c26c7759412ffcaf0dd6e3172e55d783dd7610b"
uuid = "80f14c24-f653-4e6a-9b94-39d6b0f70001"
version = "4.1.3"

[[deps.AbstractPPL]]
deps = ["AbstractMCMC", "DensityInterface", "Setfield", "SparseArrays"]
git-tree-sha1 = "6320752437e9fbf49639a410017d862ad64415a5"
uuid = "7a57a42e-76ec-4ea3-a279-07e840d6d9cf"
version = "0.5.2"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "52b3b436f8f73133d7bc3a6c71ee7ed6ab2ab754"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.3"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "195c5505521008abea5aee4f96930717958eac6f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.4.0"

[[deps.AdvancedHMC]]
deps = ["AbstractMCMC", "ArgCheck", "DocStringExtensions", "InplaceOps", "LinearAlgebra", "ProgressMeter", "Random", "Requires", "Setfield", "Statistics", "StatsBase", "StatsFuns", "UnPack"]
git-tree-sha1 = "0091e2e4d0a7125da0e3ad8c7dbff9171a921461"
uuid = "0bf59076-c3b1-5ca4-86bd-e02cd72cde3d"
version = "0.3.6"

[[deps.AdvancedMH]]
deps = ["AbstractMCMC", "Distributions", "Random", "Requires"]
git-tree-sha1 = "d7a7dabeaef34e5106cdf6c2ac956e9e3f97f666"
uuid = "5b7e9947-ddc0-4b3f-9b55-0d8042f74170"
version = "0.6.8"

[[deps.AdvancedPS]]
deps = ["AbstractMCMC", "Distributions", "Libtask", "Random", "StatsFuns"]
git-tree-sha1 = "9ff1247be1e2aa2e740e84e8c18652bd9d55df22"
uuid = "576499cb-2369-40b2-a588-c64705576edc"
version = "0.3.8"

[[deps.AdvancedVI]]
deps = ["Bijectors", "Distributions", "DistributionsAD", "DocStringExtensions", "ForwardDiff", "LinearAlgebra", "ProgressMeter", "Random", "Requires", "StatsBase", "StatsFuns", "Tracker"]
git-tree-sha1 = "67fcc7d46c26250e89fc62798fbe07b5ee264c6f"
uuid = "b5ca4192-6429-45e5-a2d9-87aec30a685c"
version = "0.1.6"

[[deps.ArgCheck]]
git-tree-sha1 = "a3a402a35a2f7e0b87828ccabbd5ebfbebe356b4"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Arpack]]
deps = ["Arpack_jll", "Libdl", "LinearAlgebra", "Logging"]
git-tree-sha1 = "9b9b347613394885fd1c8c7729bfc60528faa436"
uuid = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
version = "0.5.4"

[[deps.Arpack_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "5ba6c757e8feccf03a1554dfaf3e26b3cfc7fd5e"
uuid = "68821587-b530-5797-8361-c406ea357684"
version = "3.5.1+1"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "c46fb7dd1d8ca1d213ba25848a5ec4e47a1a1b08"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.26"

[[deps.ArrayInterfaceStaticArraysCore]]
deps = ["Adapt", "ArrayInterfaceCore", "LinearAlgebra", "StaticArraysCore"]
git-tree-sha1 = "93c8ba53d8d26e124a5a8d4ec914c3a16e6a0970"
uuid = "dd5226c6-a4d4-4bc7-8575-46859f9c95b9"
version = "0.1.3"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "1dd4d9f5beebac0c03446918741b1a03dc5e5788"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.6"

[[deps.BangBang]]
deps = ["Compat", "ConstructionBase", "Future", "InitialValues", "LinearAlgebra", "Requires", "Setfield", "Tables", "ZygoteRules"]
git-tree-sha1 = "7fe6d92c4f281cf4ca6f2fba0ce7b299742da7ca"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.3.37"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[deps.Bijectors]]
deps = ["ArgCheck", "ChainRulesCore", "ChangesOfVariables", "Compat", "Distributions", "Functors", "InverseFunctions", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "MappedArrays", "Random", "Reexport", "Requires", "Roots", "SparseArrays", "Statistics"]
git-tree-sha1 = "a3704b8e5170f9339dff4e6cb286ad49464d3646"
uuid = "76274a88-744f-5084-9051-94815aaf08c4"
version = "0.10.6"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "c5fd7cd27ac4aed0acf4b73948f0110ff2a854b2"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.7"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRules]]
deps = ["Adapt", "ChainRulesCore", "Compat", "Distributed", "GPUArraysCore", "IrrationalConstants", "LinearAlgebra", "Random", "RealDot", "SparseArrays", "Statistics", "StructArrays"]
git-tree-sha1 = "0c8c8887763f42583e1206ee35413a43c91e2623"
uuid = "082447d4-558c-5d27-93f4-14fc19e9eca2"
version = "1.45.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e7ff6cadf743c098e08fca25c91103ee4303c9bb"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.6"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "Random", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "64df3da1d2a26f4de23871cd1b6482bb68092bd5"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.14.3"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "cc4bd91eba9cdbbb4df4746124c22c0832a460d6"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.1.1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random", "SnoopPrecompile"]
git-tree-sha1 = "aa3edc8f8dea6cbfa176ee12f7c2fc82f0608ed3"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.20.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonMark]]
deps = ["Crayons", "JSON", "URIs"]
git-tree-sha1 = "86cce6fd164c26bad346cc51ca736e692c9f553c"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.8.7"

[[deps.CommonSolve]]
git-tree-sha1 = "9441451ee712d1aec22edad62db1a9af3dc8d852"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.3"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "00a2cccc7f098ff3b66806862d275ca3db9e6e5a"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.5.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.CompositionsBase]]
git-tree-sha1 = "455419f7e328a1a2493cabc6428d79e951349769"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.1"

[[deps.ConsoleProgressMonitor]]
deps = ["Logging", "ProgressMeter"]
git-tree-sha1 = "3ab7b2136722890b9af903859afcf457fa3059e8"
uuid = "88cd18e8-d9cc-4ea6-8889-5259c0d15c8b"
version = "0.1.2"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fb21ddd70a051d882a1686a5a550990bbe371a95"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.4.1"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "e08915633fcb3ea83bf9d6126292e5bc5c739922"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.13.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SnoopPrecompile", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "d4f69885afa5e6149d0cab3818491565cf41446d"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.4.4"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.DataValues]]
deps = ["DataValueInterfaces", "Dates"]
git-tree-sha1 = "d88a19299eba280a6d062e135a43f00323ae70bf"
uuid = "e7dc6d0d-1eca-5fa6-8ad6-5aecde8b7ea5"
version = "0.4.13"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DefineSingletons]]
git-tree-sha1 = "0fba8b706d0178b4dc7fd44a96a92382c9065c2c"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.2"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "c5b6685d53f933c11404a3ae9822afe30d522494"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.12.2"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "a7756d098cbabec6b3ac44f369f74915e8cfd70a"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.79"

[[deps.DistributionsAD]]
deps = ["Adapt", "ChainRules", "ChainRulesCore", "Compat", "DiffRules", "Distributions", "FillArrays", "LinearAlgebra", "NaNMath", "PDMats", "Random", "Requires", "SpecialFunctions", "StaticArrays", "StatsBase", "StatsFuns", "ZygoteRules"]
git-tree-sha1 = "0c139e48a8cea06c6ecbbec19d3ebc5dcbd7870d"
uuid = "ced4e74d-a319-5a8a-b0ac-84af2272839c"
version = "0.6.43"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "c36550cb29cbe373e95b3f40486b9a4148f89ffd"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.2"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.DynamicPPL]]
deps = ["AbstractMCMC", "AbstractPPL", "BangBang", "Bijectors", "ChainRulesCore", "ConstructionBase", "Distributions", "DocStringExtensions", "LinearAlgebra", "MacroTools", "OrderedCollections", "Random", "Setfield", "Test", "ZygoteRules"]
git-tree-sha1 = "9a795bb2fe860e2b0a19136429a173de9f8c2774"
uuid = "366bfd00-2699-11ea-058f-f148b4cae6d8"
version = "0.21.3"

[[deps.EllipticalSliceSampling]]
deps = ["AbstractMCMC", "ArrayInterfaceCore", "Distributions", "Random", "Statistics"]
git-tree-sha1 = "4cda4527e990c0cc201286e0a0bfbbce00abcfc2"
uuid = "cad2338a-1db2-11e9-3401-43bc07c9ede2"
version = "1.0.0"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.ExprTools]]
git-tree-sha1 = "56559bbef6ca5ea0c0818fa5c90320398a6fbf8d"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.8"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "90630efff0894f8142308e334473eba54c433549"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.5.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "e27c4ebe80e8699540f2d6c805cc12203b614f12"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.20"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "802bfc139833d2ba893dd9e62ba1767c88d708ae"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.5"

[[deps.FiniteDiff]]
deps = ["ArrayInterfaceCore", "LinearAlgebra", "Requires", "Setfield", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "04ed1f0029b6b3af88343e439b995141cb0d0b8d"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.17.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "187198a4ed8ccd7b5d99c41b69c679269ea2b2d4"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.32"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "a5e6e7f12607e90d71b09e6ce2c965e41b337968"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.1"

[[deps.Functors]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "a2657dd0f3e8a61dbe70fc7c122038bd33790af5"
uuid = "d9f16b24-f501-4c13-a1f2-28368ffc5196"
version = "0.3.0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "6872f5ec8fd1a38880f027a26739d42dcda6691f"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.2"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "051072ff2accc6e0e87b708ddee39b18aa04a0bc"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.71.1"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "501a4bf76fd679e7fcd678725d5072177392e756"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.71.1+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "fb83fbe02fe57f2c068013aa94bcdf6760d3a7a7"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+1"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "Dates", "IniFile", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "e1acc37ed078d99a714ed8376446f92a5535ca65"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.5.5"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions", "Test"]
git-tree-sha1 = "709d864e3ed6e3545230601f94e11ebc65994641"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.11"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "0cf92ec945125946352f3d46c96976ab972bde6f"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.3.2"

[[deps.InplaceOps]]
deps = ["LinearAlgebra", "Test"]
git-tree-sha1 = "50b41d59e7164ab6fda65e71049fee9d890731ff"
uuid = "505f98c9-085e-5b2c-8e89-488be7bf1f34"
version = "0.3.0"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "842dd89a6cb75e02e85fdd75c760cdc43f5d6863"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.6"

[[deps.IntervalSets]]
deps = ["Dates", "Random", "Statistics"]
git-tree-sha1 = "16c0cc91853084cb5f58a78bd209513900206ce6"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.4"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.InvertedIndices]]
git-tree-sha1 = "82aec7a3dd64f4d9584659dc0b62ef7db2ef3e19"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.2.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IterableTables]]
deps = ["DataValues", "IteratorInterfaceExtensions", "Requires", "TableTraits", "TableTraitsUtils"]
git-tree-sha1 = "70300b876b2cebde43ebc0df42bc8c94a144e1b4"
uuid = "1c8ee90f-4401-5389-894e-7a04a3dc0f4d"
version = "1.0.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "f377670cda23b6b7c1c0b3893e37451c5c1a2185"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.5"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "a79c4cf60cc7ddcdcc70acbb7216a5f9b4f8d188"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.16"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "9816b296736292a80b9a3200eb7fbb57aaa3917a"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.5"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LRUCache]]
git-tree-sha1 = "d862633ef6097461037a00a13f709a62ae4bdfdd"
uuid = "8ac3fa9e-de4c-5943-b1dc-09c6b5f20637"
version = "1.4.0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "ab9aa169d2160129beb241cb2750ca499b4e90e9"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.17"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LeftChildRightSiblingTrees]]
deps = ["AbstractTrees"]
git-tree-sha1 = "fb6803dafae4a5d62ea5cab204b1e657d9737e7f"
uuid = "1d6d02ad-be62-4b6b-8a6d-2f90e265016e"
version = "0.2.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtask]]
deps = ["FunctionWrappers", "LRUCache", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "dfa6c5f2d5a8918dd97c7f1a9ea0de68c2365426"
uuid = "6f1fad26-d15e-5dc8-ae53-837a1d7b8c9f"
version = "0.7.5"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "7bbea35cec17305fc70a0e5b4641477dc0789d9d"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogDensityProblems]]
deps = ["ArgCheck", "DocStringExtensions", "Random", "Requires", "UnPack"]
git-tree-sha1 = "c3e1189191e4528b605070972d7d4e9cd91dd96b"
uuid = "6fdf6af0-433a-55f7-b3ed-c6c6e0b8df7c"
version = "1.0.3"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "946607f84feb96220f480e0422d3484c49c00239"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.19"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "5d4d2d9904227b8bd66386c1138cf4d5ffa826bf"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "0.4.9"

[[deps.LombScargle]]
deps = ["FFTW", "LinearAlgebra", "Measurements", "Random", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "d64a0ce7539181136a85fd8fe4f42626387f0f26"
uuid = "fc60dff9-86e7-5f2f-a8a0-edeadbb75bd9"
version = "1.0.3"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "dedbebe234e06e1ddad435f5c6f4b85cd8ce55f7"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.2.2"

[[deps.MCMCChains]]
deps = ["AbstractMCMC", "AxisArrays", "Compat", "Dates", "Distributions", "Formatting", "IteratorInterfaceExtensions", "KernelDensity", "LinearAlgebra", "MCMCDiagnosticTools", "MLJModelInterface", "NaturalSort", "OrderedCollections", "PrettyTables", "Random", "RecipesBase", "Serialization", "Statistics", "StatsBase", "StatsFuns", "TableTraits", "Tables"]
git-tree-sha1 = "f5f347b828fd95ece7398f412c81569789361697"
uuid = "c7f686f2-ff18-58e9-bc7b-31028e88f75d"
version = "5.5.0"

[[deps.MCMCDiagnosticTools]]
deps = ["AbstractFFTs", "DataAPI", "Distributions", "LinearAlgebra", "MLJModelInterface", "Random", "SpecialFunctions", "Statistics", "StatsBase", "Tables"]
git-tree-sha1 = "dab8e9e9cd714fd3adc2344a97504e7e64e78546"
uuid = "be115224-59cd-429b-ad48-344e309966f0"
version = "0.1.5"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "2ce8695e1e699b68702c03402672a69f54b8aca9"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.2.0+0"

[[deps.MLJModelInterface]]
deps = ["Random", "ScientificTypesBase", "StatisticalTraits"]
git-tree-sha1 = "c8b7e632d6754a5e36c0d94a4b466a5ba3a30128"
uuid = "e80e1ace-859a-464e-9ed9-23947d8ae3ea"
version = "1.8.0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.MappedArrays]]
git-tree-sha1 = "e8b359ef06ec72e8c030463fe02efe5527ee5142"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.1"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MarkdownLiteral]]
deps = ["CommonMark", "HypertextLiteral"]
git-tree-sha1 = "0d3fa2dd374934b62ee16a4721fe68c418b92899"
uuid = "736d6165-7244-6769-4267-6b50796e6954"
version = "0.1.1"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Measurements]]
deps = ["Calculus", "LinearAlgebra", "Printf", "RecipesBase", "Requires"]
git-tree-sha1 = "12950d646ce04fb2e89ba5bd890205882c3592d7"
uuid = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
version = "2.8.0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.MicroCollections]]
deps = ["BangBang", "InitialValues", "Setfield"]
git-tree-sha1 = "4d5917a26ca33c66c8e5ca3247bd163624d35493"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.1.3"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.MultivariateStats]]
deps = ["Arpack", "LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI", "StatsBase"]
git-tree-sha1 = "efe9c8ecab7a6311d4b91568bd6c88897822fabe"
uuid = "6f286f6a-111f-5878-ab1e-185364afe411"
version = "0.10.0"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NNlib]]
deps = ["Adapt", "ChainRulesCore", "LinearAlgebra", "Pkg", "Requires", "Statistics"]
git-tree-sha1 = "37596c26f107f2fd93818166ed3dab1a2e6b2f05"
uuid = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"
version = "0.8.11"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

[[deps.NamedArrays]]
deps = ["Combinatorics", "DataStructures", "DelimitedFiles", "InvertedIndices", "LinearAlgebra", "Random", "Requires", "SparseArrays", "Statistics"]
git-tree-sha1 = "2fd5787125d1a93fbe30961bd841707b8a80d75b"
uuid = "86f7a689-2022-50b4-a561-43c23ac3c673"
version = "0.9.6"

[[deps.NaturalSort]]
git-tree-sha1 = "eda490d06b9f7c00752ee81cfa451efe55521e21"
uuid = "c020b1a1-e9b0-503a-9c33-f039bfc54a85"
version = "1.0.0"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "440165bf08bc500b8fe4a7be2dc83271a00c0716"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.12"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.Observables]]
git-tree-sha1 = "6862738f9796b3edc1c09d0890afce4eca9e7e93"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.4"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "f71d8950b724e9ff6110fc948dff5a329f901d64"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.8"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "df6830e37943c7aaa10023471ca47fb3065cc3c4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.3.2"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6e9dba33f9f2c44e08a020b0caf6903be540004"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.19+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "1903afc76b7d01719d9c30d3c7d501b61db96721"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.7.4"

[[deps.Optimisers]]
deps = ["ChainRulesCore", "Functors", "LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "f1cccb9f879dd4eaa4d92b115ab793545965d763"
uuid = "3bd65402-5787-11e9-1adc-39752487f4e2"
version = "0.2.13"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "cf494dca75a69712a72b80bc48f59dcf3dea63ec"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.16"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "b64719e8b4504983c7fca6cc9db3ebc8acc2a4d6"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "SnoopPrecompile", "Statistics"]
git-tree-sha1 = "5b7690dd212e026bbab1860016a6601cb077ab66"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.2"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SnoopPrecompile", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "5a554627361326403e2bb2db717ada24ae6cefbc"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.37.1"

[[deps.PlutoHooks]]
deps = ["InteractiveUtils", "Markdown", "UUIDs"]
git-tree-sha1 = "072cdf20c9b0507fdd977d7d246d90030609674b"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0774"
version = "0.0.5"

[[deps.PlutoLinks]]
deps = ["FileWatching", "InteractiveUtils", "Markdown", "PlutoHooks", "Revise", "UUIDs"]
git-tree-sha1 = "8f5fa7056e6dcfb23ac5211de38e6c03f6367794"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0420"
version = "0.1.6"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "LaTeXStrings", "Latexify", "Markdown", "PlutoLinks", "PlutoUI", "Random"]
git-tree-sha1 = "ea3e4ac2e49e3438815f8946fa7673b658e35bdb"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.2.5"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eadad7b14cf046de6eb41f13c9275e5aa2711ab6"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.49"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "LaTeXStrings", "Markdown", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "96f6db03ab535bdb901300f88335257b0018689d"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.2.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressLogging]]
deps = ["Logging", "SHA", "UUIDs"]
git-tree-sha1 = "80d919dee55b9c50e8d9e2da5eeafff3fe58b539"
uuid = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
version = "0.1.4"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "d7a7aef8f8f2d537104f170139553b14dfe39fe9"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.2"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "97aa253e65b784fd13e83774cadc95b38011d734"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.6.0"

[[deps.Query]]
deps = ["DataValues", "IterableTables", "MacroTools", "QueryOperators", "Statistics"]
git-tree-sha1 = "a66aa7ca6f5c29f0e303ccef5c8bd55067df9bbe"
uuid = "1a8c2f83-1ff3-5112-b086-8aa67b057ba1"
version = "1.0.0"

[[deps.QueryOperators]]
deps = ["DataStructures", "DataValues", "IteratorInterfaceExtensions", "TableShowUtils"]
git-tree-sha1 = "911c64c204e7ecabfd1872eb93c49b4e7c701f02"
uuid = "2aef5ad7-51ca-5a8f-8e88-e75cf067b44b"
version = "0.9.3"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[deps.RealDot]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9f0a1b71baaf7650f4fa8a1d168c7fb6ee41f0c9"
uuid = "c1ae055f-0cd5-4b69-90a6-9a35b1a98df9"
version = "0.1.0"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "18c35ed630d7229c5584b945641a73ca83fb5213"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.2"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase", "SnoopPrecompile"]
git-tree-sha1 = "e974477be88cb5e3040009f3767611bc6357846f"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.11"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterfaceCore", "ArrayInterfaceStaticArraysCore", "ChainRulesCore", "DocStringExtensions", "FillArrays", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "Tables", "ZygoteRules"]
git-tree-sha1 = "a5ce741acddc02f0d4fc6505463ca89697d7fb23"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.32.3"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "dad726963ecea2d8a81e26286f625aee09a91b7c"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.4.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.Roots]]
deps = ["ChainRulesCore", "CommonSolve", "Printf", "Setfield"]
git-tree-sha1 = "a3db467ce768343235032a1ca0830fc64158dadf"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.0.8"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "50314d2ef65fce648975a8e80ae6d8409ebbf835"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.5"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.SciMLBase]]
deps = ["ArrayInterfaceCore", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Preferences", "RecipesBase", "RecursiveArrayTools", "RuntimeGeneratedFunctions", "StaticArraysCore", "Statistics", "Tables"]
git-tree-sha1 = "6a5c8e335e82b0c674bf74f7b45f005175b0cc5f"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.77.0"

[[deps.ScientificTypesBase]]
git-tree-sha1 = "a8e18eb383b5ecf1b5e6fc237eb39255044fd92b"
uuid = "30f210dd-8aff-4c5f-94ba-8e64358c1161"
version = "3.0.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "efd23b378ea5f2db53a55ae53d3133de4e080aa9"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.16"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "38d88503f695eb0301479bc9b0d4320b378bafe5"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.8.2"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SnoopPrecompile]]
git-tree-sha1 = "f604441450a3c0569830946e5b33b78c928e1a85"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.1"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "e08a62abc517eb79667d0a29dc08a3b589516bb5"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.15"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "ffc098086f35909741f71ce21d03dadf0d2bfa76"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.11"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.StatisticalTraits]]
deps = ["ScientificTypesBase"]
git-tree-sha1 = "30b9236691858e13f167ce829490a68e1a597782"
uuid = "64bff920-2084-43da-a3e6-9bb72801c0c9"
version = "3.2.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "ab6083f09b3e617e34a956b43e9d51b824206932"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.1.1"

[[deps.StatsPlots]]
deps = ["AbstractFFTs", "Clustering", "DataStructures", "DataValues", "Distributions", "Interpolations", "KernelDensity", "LinearAlgebra", "MultivariateStats", "NaNMath", "Observables", "Plots", "RecipesBase", "RecipesPipeline", "Reexport", "StatsBase", "TableOperations", "Tables", "Widgets"]
git-tree-sha1 = "e0d5bc26226ab1b7648278169858adcfbd861780"
uuid = "f3b207a7-027a-5e70-b257-86293d7955fd"
version = "0.15.4"

[[deps.StringManipulation]]
git-tree-sha1 = "46da2434b41f41ac3594ee9816ce5541c6096123"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.0"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArraysCore", "Tables"]
git-tree-sha1 = "13237798b407150a6d2e2bce5d793d7d9576e99e"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.13"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableOperations]]
deps = ["SentinelArrays", "Tables", "Test"]
git-tree-sha1 = "e383c87cf2a1dc41fa30c093b2a19877c83e1bc1"
uuid = "ab02a1b2-a7df-11e8-156e-fb1833f50b87"
version = "1.2.0"

[[deps.TableShowUtils]]
deps = ["DataValues", "Dates", "JSON", "Markdown", "Test"]
git-tree-sha1 = "14c54e1e96431fb87f0d2f5983f090f1b9d06457"
uuid = "5e66a065-1f0a-5976-b372-e0b8c017ca10"
version = "0.2.5"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.TableTraitsUtils]]
deps = ["DataValues", "IteratorInterfaceExtensions", "Missings", "TableTraits"]
git-tree-sha1 = "78fecfe140d7abb480b53a44f3f85b6aa373c293"
uuid = "382cd787-c1b6-5bf2-a167-d5b971a19bda"
version = "1.0.2"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "c79322d36826aa2f4fd8ecfa96ddb47b174ac78d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.TerminalLoggers]]
deps = ["LeftChildRightSiblingTrees", "Logging", "Markdown", "Printf", "ProgressLogging", "UUIDs"]
git-tree-sha1 = "f53e34e784ae771eb9ccde4d72e578aa453d0554"
uuid = "5d786b92-1e48-4d6f-9151-6b4477ca9bed"
version = "0.1.6"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tracker]]
deps = ["Adapt", "DiffRules", "ForwardDiff", "Functors", "LinearAlgebra", "LogExpFunctions", "MacroTools", "NNlib", "NaNMath", "Optimisers", "Printf", "Random", "Requires", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "d963aad627fd7af56fbbfee67703c2f7bfee9dd7"
uuid = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
version = "0.2.22"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "e4bdc63f5c6d62e80eb1c0043fcc0360d5950ff7"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.10"

[[deps.Transducers]]
deps = ["Adapt", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "Setfield", "SplittablesBase", "Tables"]
git-tree-sha1 = "c42fa452a60f022e9e087823b47e5a5f8adc53d5"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.75"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.Turing]]
deps = ["AbstractMCMC", "AdvancedHMC", "AdvancedMH", "AdvancedPS", "AdvancedVI", "BangBang", "Bijectors", "DataStructures", "Distributions", "DistributionsAD", "DocStringExtensions", "DynamicPPL", "EllipticalSliceSampling", "ForwardDiff", "Libtask", "LinearAlgebra", "LogDensityProblems", "MCMCChains", "NamedArrays", "Printf", "Random", "Reexport", "Requires", "SciMLBase", "Setfield", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Tracker"]
git-tree-sha1 = "8a40377bcc4b054ebdc8f680e96cd73a4a6fe2e6"
uuid = "fce5fe82-541a-59a6-adf8-730c64b5f9a0"
version = "0.22.0"

[[deps.URIs]]
git-tree-sha1 = "ac00576f90d8a259f2c9d823e91d1de3fd44d348"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.Widgets]]
deps = ["Colors", "Dates", "Observables", "OrderedCollections"]
git-tree-sha1 = "fcdae142c1cfc7d89de2d11e08721d0f2f86c98a"
uuid = "cc8bc4a8-27d6-5769-a93b-9d913e69aa62"
version = "0.6.6"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "58443b63fb7e465a8a7210828c91c08b92132dff"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.14+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "868e669ccb12ba16eaf50cb2957ee2ff61261c56"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.29.0+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╟─f53f1d9f-8be8-4dcc-a637-f39b1eaacd63
# ╟─6d090107-d691-4a68-91ad-c7dae5df35ad
# ╟─8cf8e22c-af0f-422b-89a1-109291cd749a
# ╟─8f829f16-553a-11ed-3fd9-9f77ddd265d5
# ╟─a6eedc20-0171-48c7-938f-daee0ce4f8e9
# ╟─e71bdd5e-a50e-46b2-8249-55ad7dffb789
# ╟─4a949edb-8737-4439-a1c0-14641bd99f8c
# ╟─6f000c6d-14a8-4d6d-97f8-d01f4dd516bc
# ╟─e4084b5c-19de-460d-aceb-ee6f3cfafbd9
# ╟─e3b46289-f602-4dee-81cc-91c8466dcd5a
# ╟─64987161-afde-4b4d-93a8-19f7b40aeae7
# ╟─383d9191-fc82-4f9c-81f5-837a67c71e9b
# ╟─10290d00-c667-46f3-a03e-8fb33842448a
# ╟─0b49db68-973a-42e0-8eed-4f4405b5f808
# ╟─5a929ea7-93cc-4c4d-9340-734bcd509719
# ╟─b5fc8991-ae2d-4db4-a7f9-5da11a21e094
# ╟─7f8751c3-b087-401e-a098-25460bf496d9
# ╟─7e9f805c-80c4-472f-a3ea-f73a732b8a57
# ╟─62e93800-22bd-40a1-9ee6-bcd1134539ae
# ╟─d63c6ac8-3623-426e-b25f-668ffce47ddd
# ╟─44935b70-11a6-4804-a7cb-150d8660441b
# ╟─09020d68-05d3-4908-9dcc-c1e0097cd54c
# ╟─e53d6766-4fb1-4bda-822e-9e468b440fdf
# ╟─1ce61116-86e3-4488-b2ec-a0d7617bb4b3
# ╟─1121402e-69d9-4d3b-bec9-c4011e226675
# ╟─b51622f7-baf7-4ab5-861e-dd7f5a707fdb
# ╟─ec77269a-6cf8-420f-a7ef-75f15e30de28
# ╟─68a49dcb-3cd9-4905-abf6-d43083d256ce
# ╟─3b21b28a-52f0-4da0-beba-bb09d666b9c6
# ╟─f46cbed8-dec6-4bb2-8942-4fe96ebea3e4
# ╟─aae7b265-1ff7-4d5b-8d93-f0bb0bce575a
# ╟─2384f27b-05d7-45f8-a18c-75bcdea14b80
# ╟─58b75f36-b9f5-4192-a796-73f65e497df2
# ╟─6b5c4ec5-3b54-47e6-ad0b-d9a813b7bb7a
# ╟─db3dca2a-86ed-4c6e-9477-175092eb538f
# ╟─69f83e46-8788-42ae-945d-3fb27b78a36f
# ╟─1d9f11ce-1515-4d75-9a30-55db3b3fa8ae
# ╟─255ce271-2011-4082-acaa-bfd4d972dc7b
# ╟─3a51761b-8659-4c42-9612-28d6daaf7404
# ╟─e8222797-77df-40d2-b9bd-279ed3c12cde
# ╟─06eede8c-1dda-4a15-a69c-6019ee11269d
# ╟─a494b98a-08df-464c-b3f4-584463f4210c
# ╟─fb7ddbe7-6f42-4a6d-8d28-611df2cdf549
# ╟─7f720fc1-3004-41af-a570-46183ed4b646
# ╟─c344a8d1-da07-4206-87d2-92924e133001
# ╟─6926b12f-d5e9-475a-8248-2a4f0e0f18b5
# ╟─22ebdf59-a508-44f3-9db8-031b37c4446d
# ╟─2d40a7b2-f650-45e0-978e-eca24954faa6
# ╟─acf3b224-3aa7-4419-92a1-fd2b1ee582ce
# ╟─a40a805d-48a7-4451-b958-56ef041d3333
# ╟─200c1f28-0eb3-459f-91f8-a0c38cb92bed
# ╟─59069f4e-58fa-458a-8c47-4176e3970f43
# ╟─17b1760d-ed77-4773-8147-43245c70a962
# ╟─5712b24d-cfff-4c95-8008-dfad4dc32a2c
# ╟─34304691-cdd8-4cc0-a33e-e166c434eb4d
# ╟─79e5a375-3c16-4443-86b9-d3bbad6101d4
# ╟─34fa0f41-fae2-4854-8055-d9ee476c3eef
# ╟─36e7d4ab-735e-4517-bf2e-ae1db69d227e
# ╟─519b083c-aa8d-48da-b9ec-6e3cddf94d99
# ╟─f3bbb76c-14cd-4566-935b-5c1f6949eaf5
# ╟─2350b4a6-a538-430e-b832-4ecb4a458c4d
# ╟─2de85bb1-881f-483a-bcdb-2109ed795ec5
# ╟─bfc1affc-b392-4884-a35f-55593af7db53
# ╟─8f239f05-2610-42ce-92c6-968943418328
# ╟─32a6f1f8-c213-40c1-a6c1-5dc7ce87976e
# ╟─a435c976-de91-40a6-b536-9cf005027a3b
# ╟─640c9851-73ba-47ed-84e4-6484f3b793b2
# ╟─59bc7b62-87f8-42bb-80d0-74d4cf2c9edf
# ╟─963230d5-cef8-4a9c-aef4-319cbd75c207
# ╟─d24fbe45-dc19-4946-9958-b3f26650d572
# ╟─75e5df52-18ab-4046-a2a0-fdaea68d9543
# ╟─3e653573-04db-4913-8771-2c23fe0e01a1
# ╟─c2c32163-dcbc-4bfc-a172-e3d750ce6c4b
# ╟─49a9a26e-6423-47c3-91ca-785c2ccafe24
# ╟─70e8d18a-fa92-4993-8994-c9a481141e1d
# ╟─acff1e9d-038c-4267-9f85-37e21122988c
# ╟─4a3292ab-cfa2-4b58-ba66-19a23f8f2a23
# ╟─9e5a3e25-8f7a-469e-aa54-846beec8990f
# ╟─7d7fba67-c78a-44fe-bd9b-d8d4967b32c7
# ╟─e30755c3-e701-4fd6-a244-6755ed6603d3
# ╟─769de9c4-6167-4550-be08-320c7c63fe3e
# ╠═7c2ceba1-6e8a-4a5f-845d-73b97810c099
# ╟─9d1eadf9-7fe9-4620-b866-4a442033a0b4
# ╟─5736df1b-a0c9-47b2-9cfe-faf01641ce26
# ╟─74170a94-2014-45fd-9f1d-d476b1febf14
# ╟─0f2ee09b-1d31-4a34-9b0a-22e8b91f4303
# ╟─6cea0d3e-daed-4318-83b0-13bf9fe00d2a
# ╟─5dcd30d7-dd60-4b15-96a6-4b222e55d779
# ╟─c67d1588-36e4-44aa-a98b-8308bf57e1e0
# ╟─d892f267-8a1f-41e8-9104-573ec424abc9
# ╟─2293d249-407a-4578-8eb8-377a3ef44e68
# ╟─71e66a2a-be4d-48f3-a868-bcd8f7647e22
# ╟─0372f0ee-b238-403e-b045-1ade617a271d
# ╟─61379c75-6969-4652-b6ce-7e1e992f1f23
# ╟─1d397e31-740e-41e7-ab4c-d6009752ee33
# ╟─514730e4-f589-4841-a4db-4b39e746e6a9
# ╟─96851df0-0aad-44ac-830a-b5062917f4cf
# ╟─167cc8ad-cb04-47a1-bb25-08e9c345b24e
# ╟─5edca12e-e708-43d3-b0dc-4b70c1c4ea70
# ╠═66446e92-0d84-4d2b-9804-a563fbd8e30b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
