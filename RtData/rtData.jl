# load packages
println("loading packages ...")
using StatsBase, Distributions, DataFrames, CSV, Dates, HTTP, Impute, DSP, ExcelReaders;
# 
# using Plots, 

# set the output url for saving data and plots
# directories called data and plots must exist at this location
output_url = "./"
clean_data = DataFrame(parent=String[], location=String[], date=Date[], cases=Float64[])
allowmissing!(clean_data)


#load global data from OWID
println("loading World data at Country level ...")
url="https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/ecdc/full_data.csv"
data = DataFrame(CSV.read(IOBuffer(HTTP.get(url).body)))
select!(data, [:location, :date, :total_cases])
insertcols!(data, 1, :parent => repeat(["World"],size(data)[1]))
rename!(data,names(clean_data))
sort!(data,(:parent,:location,:date))
clean_data = vcat(clean_data, data)
println("World data loaded, size = ",size(clean_data)) 

# load all UK data from Tom White github repository
println("loading United Kingdom data at Upper tier level ...")
url="https://raw.githubusercontent.com/tomwhite/covid-19-uk-data/master/data/covid-19-cases-uk.csv"
data = DataFrame(CSV.read(IOBuffer(HTTP.get(url).body)))
data = data[.!(ismissing.(data[:,:AreaCode])),:]
select!(data, Not(:AreaCode))
data=data[:,[:Country, :Area, :Date, :TotalCases]]
rename!(data,names(clean_data))
sort!(data,(:parent,:location,:date))
# append to clean_data
clean_data = vcat(clean_data, data)

# now do England totals
println("loading totals for England ... ")
url="https://raw.githubusercontent.com/tomwhite/covid-19-uk-data/master/data/covid-19-totals-england.csv"
data = DataFrame(CSV.read(IOBuffer(HTTP.get(url).body)))
data = data[.!(ismissing.(data[:,:ConfirmedCases])),:]
select!(data, [:Date, :ConfirmedCases])
insertcols!(data, 1, :parent => repeat(["United Kingdom"],size(data)[1]))
insertcols!(data, 2, :location => repeat(["England"],size(data)[1]))
rename!(data,names(clean_data))
sort!(data,(:parent,:location,:date))
# append to clean_data
clean_data = vcat(clean_data, data)

# now do Scotland totals
println("loading totals for Scotland ... ")
url="https://raw.githubusercontent.com/tomwhite/covid-19-uk-data/master/data/covid-19-totals-scotland.csv"
data = DataFrame(CSV.read(IOBuffer(HTTP.get(url).body)))
data = data[.!(ismissing.(data[:,:ConfirmedCases])),:]
select!(data, [:Date, :ConfirmedCases])
insertcols!(data, 1, :parent => repeat(["United Kingdom"],size(data)[1]))
insertcols!(data, 2, :location => repeat(["Scotland"],size(data)[1]))
rename!(data,names(clean_data))
sort!(data,(:parent,:location,:date))
# append to clean_data
clean_data = vcat(clean_data, data)

# now do Wales totals
println("loading totals for Wales ... ")
url="https://raw.githubusercontent.com/tomwhite/covid-19-uk-data/master/data/covid-19-totals-wales.csv"
data = DataFrame(CSV.read(IOBuffer(HTTP.get(url).body)))
data = data[.!(ismissing.(data[:,:ConfirmedCases])),:]
select!(data, [:Date, :ConfirmedCases])
insertcols!(data, 1, :parent => repeat(["United Kingdom"],size(data)[1]))
insertcols!(data, 2, :location => repeat(["Wales"],size(data)[1]))
rename!(data,names(clean_data))
sort!(data,(:parent,:location,:date))
# append to clean_data
clean_data = vcat(clean_data, data)

# now do Northern Ireland totals
println("loading totals for Northern Ireland ... ")
url="https://raw.githubusercontent.com/tomwhite/covid-19-uk-data/master/data/covid-19-totals-northern-ireland.csv"
data = DataFrame(CSV.read(IOBuffer(HTTP.get(url).body)))
data = data[.!(ismissing.(data[:,:ConfirmedCases])),:]
select!(data, [:Date, :ConfirmedCases])
insertcols!(data, 1, :parent => repeat(["United Kingdom"],size(data)[1]))
insertcols!(data, 2, :location => repeat(["Northern Ireland"],size(data)[1]))
rename!(data,names(clean_data))
sort!(data,(:parent,:location,:date))
# append to clean_data
clean_data = vcat(clean_data, data)

println("All United Kingdom data loaded, size = ",size(clean_data)) 

# now USA data from covidtracking.com
println("loading United States data at state level ...")
url="https://covidtracking.com/api/v1/states/daily.csv"
data = DataFrame(CSV.read(IOBuffer(HTTP.get(url).body)))
# dates are integers, format as dates
data[:,:format_date] = (s->Date(string(s),DateFormat("yyyymmdd"))).(data[:,:date])
select!(data, [:state, :format_date, :total])
insertcols!(data, 1, :parent => repeat(["United States"],size(data)[1]))
rename!(data,names(clean_data))
sort!(data,(:parent,:location,:date))
clean_data = vcat(clean_data, data)
println("United States data loaded, size = ",size(clean_data)) 


# load South African data from the dsfsi
println("loading South Africa data at provincial level ...")
url="https://raw.githubusercontent.com/dsfsi/covid19za/master/data/covid19za_provincial_cumulative_timeline_confirmed.csv"
data = DataFrame(CSV.read(IOBuffer(HTTP.get(url).body)))
data[:,:format_date] .= (s->Date(string(s),DateFormat("dd-mm-yyyy"))).(data[:,:date])
function extract_za_prov(data)
    data_long = DataFrame(parent=String[], location=String[], date=Date[], cases=Float64[]) 
    for prov in names(data)[3:11]
        data_prov = data[:,[:format_date, prov]]
        insertcols!(data_prov, 1, :parent => repeat(["South Africa"],size(data_prov)[1]))
        insertcols!(data_prov, 2, :location => repeat([string(prov)],size(data_prov)[1]))
        rename!(data_prov,names(data_long))
        data_long = vcat(data_long, data_prov)
        data_long = data_long[completecases(data_long), :]
    end
    return data_long 
end
clean_data = vcat(clean_data, extract_za_prov(data))
println("South Africa loaded, size = ",size(clean_data)) 

# load Italy data from pcm_dpc 
println("loading Italy data at regional level ...")
url = "https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-regioni/dpc-covid19-ita-regioni.csv"
data = DataFrame(CSV.read(IOBuffer(HTTP.get(url).body)))
# dates are integers, format as dates
data[:,:format_date] .= (s->Date(string(s)[1:10])).(data[:,:data])
select!(data, [:denominazione_regione, :format_date, :totale_casi])
insertcols!(data, 1, :parent => repeat(["Italy"],size(data)[1]))
rename!(data,names(clean_data))
data = data[completecases(data), :]
sort!(data,(:parent,:location,:date))
clean_data = vcat(clean_data, data)
println("Italy data loaded, size = ",size(clean_data)) 

# load Sweden data, xlsx format
println("loading Sweden data at regional level ...")
url="https://www.arcgis.com/sharing/rest/content/items/b5e7488e117749c19881cce45db13f7e/data"
data = DataFrame(readxlsheet(download(url), "Antal per dag region", skipstartrows=0))
data_names=Symbol.(Array(data[1,:]))
data = DataFrame(readxlsheet(download(url), "Antal per dag region", skipstartrows=1))
rename!(data,data_names)
data[:,:format_date] .= (s->Date(string(s)[1:10])).(data[:,:Statistikdatum])
ddf = DataFrame(parent=String[], location=String[], date=Date[], cases=Int64[])
for location in data_names[3:end]
    df = data[:,[:format_date, location]]
    insertcols!(df, 1, :parent => repeat(["Sweden"],size(data)[1]))
    insertcols!(df, 2, :location => repeat([string(location)],size(data)[1]))
    rename!(df,Symbol.(["parent", "location", "date", "cases"]))
    df[:,:cases]=cumsum(Float64.(df[:,:cases]))
    for j in 1:size(df)[1]
        push!(ddf, df[j,:])
    end
end
clean_data = vcat(clean_data, ddf)
println("Sweden data loaded, size = ",size(clean_data)) 


sort!(clean_data,(:parent,:location,:date))
println("All data loaded, latest data at ", maximum(clean_data[:,:date]))
println("computing Rt estimates ... ")



# some model parameters that will be used throughout this script
IFP = 7   # infectious period in days
GAMMA = 1/IFP  # recovery rate
SMOOTHING_PASSES = 3 # number of times to apply moving average filter to case data
INF_CUTOFF = 0.0 # to detect when epidemic is over
ONSET_DELAY = 4  # default onset delayin days, might get changed later
FULL_BAYES = true # set true for Bayesean computation of mlrt plus density interval
println("IFP=",IFP," GAMMA=",GAMMA," SMOOTHING_PASSES=",
    SMOOTHING_PASSES," INF_CUTOFF=",INF_CUTOFF," ONSET_DELAY=",ONSET_DELAY, 
    " FULL_BAYES=",FULL_BAYES)

# create Weibull filter and helper ftns for adjusting for reporting delay
adj_range = collect(1:1:30)
p_delay = [pdf.(Weibull(2, 5),t) for t in adj_range]
p_delay = p_delay ./ sum(p_delay)

function estimate_onset(new_cases, p_delay)
    return reverse(conv(reverse(new_cases),p_delay))[length(p_delay):end]
end

function adjust_onset(onset, p_delay)
    c_delay = cumsum(p_delay)
    extras = length(onset) - length(c_delay)
    if (extras > 0)
        cd = vcat(c_delay,ones(extras))
        adjusted = onset ./ reverse(cd)
    else
        adjusted = onset
    end
    return adjusted[1:end]
end

function moving_average(inf; alpha=0.5)
    vs = inf
    if length(vs) > 1
        ret = 
            # vcat( vs[1], [ alpha * vs[i] + (1 - alpha)*vs[i-1] for i = 2:length(vs) ] )
            vcat((vs[1]+2*vs[1])/3,
                   [sum(@view vs[i:(i+2)])/3 for i in 1:(length(vs)-(2))],
                    (2*vs[length(vs)]+vs[length(vs)])/3)
    else
        ret = vs
    end
    return ret
    return Int64.(round.(ret))
end

# function to compute size of infectious pool from adjusted case counts
function compute_infectives(new_cases) 
    onset = estimate_onset(new_cases,p_delay)
    adj_onset = adjust_onset(onset,p_delay)
    if length(adj_onset) >= IFP
       infects = [sum(adj_onset[i-IFP+1:i])
                        for i in IFP:length(adj_onset) ] 
    else
        infects = []
    end
    return infects
end

# function to process data from a single location (ie one time series)
function process_data(data)
    nz = copy(data[data[:,:date] .> Date(2020,01,01),:])
    sort!(nz,:date)
    allowmissing!(nz)
    # now insert missing dates and impute cumulative cases for those dates
    nz = nz[completecases(nz), :]
    if size(nz)[1] > 1
        dr = nz[1,:date]:Day(1):nz[end,:date]
        for d in dr
            if ! (d in nz[:,:date])
                push!(nz,[missing for i in 1:size(nz)[2]])
                nz[end,:location] = nz[end-1,:location]
                nz[end,:parent] = nz[end-1,:parent]
                nz[end,:date] = d
            end
        end
        sort!(nz,:date)
        nz[!,:cases] = Impute.interp(nz[:,:cases])
    end
    smooth_new_cases = diff(nz[:,:cases])
    for i in 1:SMOOTHING_PASSES
        smooth_new_cases = moving_average(smooth_new_cases)
    end
    infectives = compute_infectives(smooth_new_cases)
    nz = nz[IFP+1:end,:]
    nz[:,:infectives] = Float64.(infectives)
    return nz   
end 

function strip_leading_zero_case_counts(data)
    nz = copy(data)
    nc = nz[:,:cases]
    ind = 1
    while (ind <= length(nc) && ( nc[ind]==0 || isnan(nc[ind])) )
        ind+=1
    end
    if (ind != 1)
        nz = copy(nz[ind:end,:])
    end
    # make sure cases are increasing
    ind = 2
    while (ind <= size(nz)[1])
        if nz[ind,:cases] < nz[ind-1,:cases]
            nz[ind,:cases] = nz[ind-1,:cases]
        end
        ind+=1
    end
    return(nz)
end

# function to extract data from one location, strip zeros and process
function prepare_cases(full_data, parent, location)
    location_data = full_data[full_data[:,:parent].==parent,:]
    location_data = location_data[location_data[:,:location].==location,:]
    sort!(location_data,:date) 
    location_data = strip_leading_zero_case_counts(location_data)
    if size(location_data)[1] < 10  # need at least one month of case counts
        return (location_data, false)
    end
    location_data_filtered = process_data(location_data)
    sanity_check=true
    if (ismissing(minimum(location_data_filtered[:,:infectives])) ||
        minimum(location_data_filtered[:,:infectives])< 0)
        sanity_check = false
    end
    return (location_data_filtered, sanity_check)
end

# function to estimate Rt from It for data from a single location that has already been pre-processed.
function estimate_rt(data)
    infectives = data[:,:infectives]
    # instantaneous Rt estimate
    rt_est = [ 1 + IFP * ((infectives[i+1]-infectives[i])/(infectives[i])) 
                            for i in 1:length(infectives)-1]
    # check for divide by zero
    rt_est .= ifelse.(isnan.(rt_est), 0, rt_est)
    # apply low Rt cutoff 
    for i in 1:length(rt_est)
        max_inf = maximum(infectives[1:i])
        if (( infectives[i]/max_inf < INF_CUTOFF ) || rt_est[i]<0 )
            rt_est[i] = 0
        end
    end
    # sanity check
    sane = true
    last_week = rt_est[end-6:end]
    if maximum(last_week)-minimum(last_week) > 100
        sane = false
        print("something fishy...")
    end
    data = data[1:end-1,:]
    data[!,:rt_est] = rt_est
    # throw away the last ONSET_DELAY data
    data = data[1:end-ONSET_DELAY+1,:]
    return (data,sane)
end     

#
# now for a whole lot of functions for Bayesean inference 
#

function highest_density_interval(pmf, p=.9)
    fb = findmax(pmf,dims=1)[2][1]
    cs = cumsum(pmf, dims=1)
    first = [1, length(cs)]
    best = first
    for (i, value) in enumerate(cs)
        for (j, high_value) in enumerate(cs[i+1:end])
            if ( (high_value-value > p) && (j<best[2]-best[1]) )
                best = [i, i+j]
                break
            end
        end
    end
    if (best == first) 
        # println(" invoking fallback ... ", fb )
        best = [fb, fb]
    end
    return best
end    

function get_posteriors(inf, sigma=0.23)
    sr = inf # Int64.(floor.(inf))
    # (1) Calculate Lambda
    RT_MAX = 12
    rt_range = collect(range(0, RT_MAX, length=RT_MAX*100+1))
    # GAMMA = 1/IFP   # done previously
    lambdas = [ sr[i-1] .* exp.(GAMMA .* (rt_range .- 1)) for i in 2:length(sr)]
    lambdas = hcat(lambdas...)

    # (2) Calculate each day's likelihood
    pDist = [[Poisson(lam) for lam in lambdas[:,j]] for j in 1:size(lambdas)[2] ]
    lks = [ [pdf.(pDist[j-1][i],Int64(floor(sr[j]))) for i in 1:length(pDist[j-1])] for j in 2:length(sr)]
    lks = hcat(lks...)
    likelihoods = lks ./ sum(lks, dims=1)
    
     # (3) Create the Gaussian Matrix
    nDist = [ Normal(rt,sigma) for rt in rt_range ]
    process_matrix = [ [pdf.(nDist[j],rt_range[i]) for i in 1:length(rt_range)] 
        for j in 1:length(rt_range)]
    process_matrix = hcat(process_matrix...)
    
    # (3a) Normalize all rows to sum to 1
    process_matrix = process_matrix ./ sum(process_matrix, dims=2)
    
    # process_matrix = I   # if you want current prior = previous posterior
    
    # (4) Calculate the initial prior
    guess0 = 2+log(sr[2]/sr[1])/GAMMA
    if guess0 < 1 
        guess0 = 1
    end
    # println(guess0)
    prior0 = [ pdf.(Gamma(guess0),rt) for rt in rt_range ]
    prior0 = prior0 ./ sum(prior0)
    
    # Create a DataFrame that will hold our posteriors for each day
    # Insert our prior as the first posterior.
    posteriors = DataFrame( [ Real for i in 1:length(sr) ], 
        [ Symbol.("D$i") for i in 1:length(sr) ], length(rt_range) )
    days = names(posteriors)
    # posteriors[:,:rts] = rt_range
    posteriors[:,:D1] = prior0
    
    # We said we'd keep track of the sum of the log of the probability
    # of the data for maximum likelihood calculation.
    log_likelihood = 0.0
    
       # (5) Iteratively apply Bayes' rule
    nd = length(days)
    # println("nd=$nd")
    for (previous_day, current_day) in zip(1:nd-1, 2:nd)

        #(5a) Calculate the new prior
        current_prior =   process_matrix * posteriors[:,Symbol("D$previous_day")] # 
        
        #(5b) Calculate the numerator of Bayes' Rule: P(k|R_t)P(R_t)
        numerator = likelihoods[:,current_day-1] .* current_prior
        
        #(5c) Calcluate the denominator of Bayes' Rule P(k)
        denominator = sum(numerator)
        
        # Execute full Bayes' Rule
        posteriors[:,Symbol("D$current_day")] = numerator./denominator
        
        # Add to the running sum of log likelihoods
        log_likelihood += log(denominator)
        
    end
    
    return (posteriors, rt_range, log_likelihood)
    
end

# the rest of the script loops through each location in the data set
# calculates Rt and It values at each date, gathers the data and at the end creates a static
# ranking plot and outputs the data as a csv for dynamic display

parents = unique(clean_data[:,:parent])
ddf = DataFrame()
if FULL_BAYES
    ddf = DataFrame(parent=String[], location=String[], date=Date[], 
            ct=Int64[], it=Int64[], rt_est=Float64[], rt_ml=Float64[], low=Float64[], high=Float64[])
else   
    ddf = DataFrame(parent=String[], location=String[], date=Date[], ct=Int64[], it=Int64[], rt_est=Float64[])
end
allowmissing!(ddf)
for parent in parents 
    cdf = clean_data[clean_data[:,:parent].==parent,:]
    locations = unique(cdf[:,:location])
    for location in locations
        dc,sane = prepare_cases(clean_data,parent,location)
        print(parent," -> ")
        if ( ! sane ) 
            println("skipping $location ... sanity check is $sane")
        else
            infectives = dc[:,:infectives]
            if ( length(infectives)<10 || maximum(infectives) < 20)
                println("skipping $location ... not enough data")
            else
                println("processing $location ...")
                # now use helpter to compute rt estimates
                dc,sane = estimate_rt(dc)
                if (! sane)
                    println("skipping $location ... Rt estimation failed")
                else
                    est_date = dc[end,:date]
                    # do full Bayesean inference if flag set
                    if FULL_BAYES
                        (posteriors, rts, log_likelihood) = get_posteriors(dc[:,:infectives])
                        mlrt_values = rts[map(x->x[1],findmax(Array(posteriors),dims=1)[2])[1,:]]
                        inds = [ highest_density_interval(posteriors[:,k]) for k in 1:size(posteriors)[2] ] 
                        low_high = [rts[inds[k]] for k in 1:length(inds)]
                        lows = map(x->x[1],low_high)
                        highs = map(x->x[2],low_high)
                        dc[:,:rt_ml]=mlrt_values
                        dc[:,:low]=lows
                        dc[:,:high]=highs
                    end
                    # now collect all the estimated data into one dataframe
                    for dp in 1:length(dc[:,:date])
                        if FULL_BAYES
                            push!(ddf,[dc[dp, :parent] location dc[dp,:date] Int64(round(dc[dp,:cases])) Int64(round(dc[dp,:infectives])) round(dc[dp,:rt_est],digits=2) round(dc[dp,:rt_ml],digits=2) round(dc[dp,:low],digits=2) round(dc[dp,:high],digits=2) ])
                        else
                            push!(ddf,[dc[dp, :parent] location dc[dp,:date] Int64(round(dc[dp,:cases])) Int64(round(dc[dp,:infectives])) round(dc[dp,:rt_est],digits=2) ])
                        end
                    end
                end
            end
        end
    end
end

# save all the estimates as one big csv file for charting
sort!(ddf,(:parent,:location,:date))
CSV.write(output_url*"data/rt_estimates.csv",ddf)

println("All done, latest Rt estimate at ", maximum(ddf[:,:date]))

println("Sample plot ... ")

# create a sample plot 
using Plots
location="South Africa"
df = ddf[ddf[:,:location].==location,:]
est_date = maximum(df[:,:date])
est_score = df[end,:rt_est]
fig_rts = plot(legendtitle="latest R from:",legend=:topright,ylabel="R(t)", ylimits = [-1,4],
    title="tracking R(t) in $location \n up to $est_date")
plot!(fig_rts, df[:,:date], df[:,:rt_est], lw=6, color=1,label="Rt_est $est_score")
if FULL_BAYES
    mlrt_score = df[end,:rt_ml]
    plot!(fig_rts, df[:,:date], df[:,:rt_ml], lw=6, color=3,
        ribbon = (df[:,:rt_ml]-df[:,:low],df[:,:high]-df[:,:rt_ml]),
        fillalpha=0.2,label="Rt_ml $mlrt_score")
end
plot!(fig_rts, df[:,:date], ones(size(df)[1]), lw=6, color=2, label="Rt=1")
savefig(fig_rts,output_url*"plots/RtLiveSample")

println("sample plot completed")

exit()


#=
# first attempt at getting data from the 4 UK nations before I found Tom White github
#
# load Scotland data from XLSX file 
println("loading Scotland data ...")
url="https://www.gov.scot/binaries/content/documents/govscot/publications/statistics/2020/04/coronavirus-covid-19-trends-in-daily-data/documents/covid-19-data-by-nhs-board/covid-19-data-by-nhs-board/govscot%3Adocument/COVID-19%2Bdata%2Bby%2BNHS%2BBoard%2B27%2BMay%2B2020.xlsx"
data_names = readxlsheet(download(url), "Table 1 - Cumulative cases", skipstartrows=2, nrows=1)[1,:]
data = DataFrame(readxlsheet(download(url), "Table 1 - Cumulative cases", skipstartrows=3))
rename!(data,Symbol.(data_names))
# loop through each cell and update clean_data
for i in 1:size(data)[1]
    for j in 2:size(data)[2]
        date = Date(data[i,1])
        parent = "Scotland"
        location = String(names(data)[j])
        if location==parent
            parent="United Kingdom"
        end
        cst = data[i,j]
        if cst != "*"
            if typeof(cst)!=Float64
                cases = parse(Float64, cst)
            else
                cases = cst
            end
            push!(clean_data,[parent location date cases])
        end
    end
end
sort!(clean_data,[:parent,:location,:date])
println("Scotland data loaded, latest at ",maximum(clean_data[:,:date]))

#
# load England data
# and pick out Upper tier plus National data
println("loading England data ...")
url="https://coronavirus.data.gov.uk/downloads/csv/coronavirus-cases_latest.csv"
full_data = CSV.read(IOBuffer(HTTP.get(url).body))
col_names=Symbol.(["Area name", "Area type", "Specimen date", "Daily lab-confirmed cases" ])
data=DataFrame(full_data[:,col_names])
rename!(data,Symbol.(["location","type","date","cases"]))
# two parts to get, first all the Upper tier wards then England at nation level
data1 = data[data[:,:type].=="Upper tier local authority",:]
parent = []
for i in 1:size(data1)[1] push!(parent,"England") end
select!(data1, Not(:type))
data1 = hcat(DataFrame(parent=parent), data1)
# now nation data 
data2 = data[data[:,:type].=="Nation",:]
parent = []
for i in 1:size(data2)[1] push!(parent,"United Kingdom") end
select!(data2, Not(:type))
data2 = hcat(DataFrame(parent=parent), data2)
# concatonate them
data = vcat(data1, data2)
data[:,:cases] = round.(data[:,:cases],digits=0)

sort!(data,[:parent,:location,:date])
locations = sort(unique(data[:,:location]))

for location in locations
    indx=data[:,:location].==location
    new_cases=data[indx,:cases]
    cases = cumsum(new_cases)
    data[indx,:cases]=cases
end
println("England data loaded, latest at ",maximum(data[:,:date]))

clean_data = vcat(clean_data, data)
sort!(clean_data,[:parent,:location,:date])
println("latest clean data at ",maximum(clean_data[:,:date]))


printstyled(first(clean_data,60))
printstyled(last(clean_data,60))

locations = sort(unique(clean_data[:,:location]))
printstyled(locations)

exit()
=#