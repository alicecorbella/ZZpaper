# Direcotries
CODEwd = "/Users/alice/projects/ZZpaper/"
DATAwd = "/Users/alice/projects/ZZpaper/simulacrum_release_v1.2.0.2017/data/"

using CSV
using DataFrames
using StatsBase
using Serialization


patients = DataFrame(CSV.File(string(DATAwd, "sim_av_patient.csv")))
tumours = DataFrame(CSV.File(string(DATAwd, "sim_av_tumour.csv")))
sact_patients = DataFrame(CSV.File(string(DATAwd, "sim_sact_patient.csv")))
sact_tumours = DataFrame(CSV.File(string(DATAwd, "sim_sact_tumour.csv")))
sact_cycle = DataFrame(CSV.File(string(DATAwd, "sim_sact_cycle.csv")))
sact_regimen = DataFrame(CSV.File(string(DATAwd, "sim_sact_regimen.csv")))
sact_outcome = DataFrame(CSV.File(string(DATAwd, "sim_sact_outcome.csv")))
sact_drugs = DataFrame(CSV.File(string(DATAwd, "sim_sact_drug_detail.csv")))
p_summ = describe(patients)
t_summ = describe(tumours)
sp_summ = describe(sact_patients)
st_summ = describe(sact_tumours)
sc_summ = describe(sact_cycle)
sr_summ = describe(sact_regimen)
so_summ = describe(sact_outcome)
sd_summ = describe(sact_drugs)



print(p_summ)
print(t_summ)
print(sp_summ)
print(st_summ)
print(sc_summ)
print(sr_summ)
print(so_summ)
print(sd_summ)



# join the two datasets on patients id
pt = leftjoin(patients, tumours, on = :PATIENTID, makeunique=true)
# # CHECK if the join is the same left and right
# pt2 = rightjoin(patients, tumours, on = :PATIENTID, makeunique=true)
# isequal(pt,pt2)
# # CHECKED 👍

pt_summ = describe(pt)
print(pt_summ)

# there are patients with multiple cancer
FT = combine(groupby(pt, :PATIENTID), nrow => :Freq)
sort(FT, [:Freq])
sum(FT[findall(FT.Freq.>1),:Freq]) - length(findall(FT.Freq.>1))
pt[findall(pt.PATIENTID.==110010131), vcat(1:3,12:18)]
# the tumour_id was given increasing order by diagnosis date
# or else if multiple cancer were diagnosed on the same day.

# get the patient id of all the people who have more than a tumour (mt1t)
pid_mt1t = unique(pt[nonunique(pt[:,[:PATIENTID]]), :PATIENTID])

# we could check the record of the first patient who had mt1t
findall(pt.PATIENTID.==pid_mt1t[1])
print(pt[findall(pt.PATIENTID.==pid_mt1t[1]), vcat(1:3,12:18)])

# create a vecotor with all the indexes to be delated
tbd_idx = Vector(undef, size(tumours)[1]-size(patients)[1])
# fill it in
j=1
for idx_mt1t_i in 1:length(pid_mt1t)
   # step 1 : select the indexes of a perso who has more than a tumor
   idxs_mt1t = findall(pt.PATIENTID.==pid_mt1t[idx_mt1t_i])
   # step 2 : select the tumorid of this patient except the smallest (first diagnosed)
   tid_tbd = sort(pt[idxs_mt1t, :TUMOURID])[2:length(idxs_mt1t)]
   # step 3 : select the lines with toumorid equal to the 2nd third tumor ...
   for nt in 1:length(tid_tbd)
      tbd_idx[j] = findall(pt.TUMOURID.==tid_tbd[nt])[1]
      j=j+1
   end
   print(string(idx_mt1t_i/length(pid_mt1t), "\n"))
end

# build the dataset of all the people with their first diagnosed cancer
pfdt = pt[Not(tbd_idx), :]

# CHECK to be sure
freq_pfdy = combine(groupby(pfdt, :PATIENTID), nrow => :Freq)
# CHECKED : all have frequencies 1 👍

# get the summaries
pfdt_summ=describe(pfdt)
print(pfdt_summ)

# WORKING ON COVARIATES
# is sex the same?
isequal(pfdt.SEX, pfdt.SEX_1)
# yes, let's drop it!
select!(pfdt, Not(:SEX_1))

# add a column if had more than one toumor
MULTIPLETUMOUR = zeros(size(pfdt)[1])
# patient id of people with multiple tumor
pid_mt = FT[findall(FT.Freq.>1),:PATIENTID]
# get the column of people with multiple tumor
MULTIPLETUMOUR = in(pid_mt).(pfdt[:,:PATIENTID])
# add the column
pfdt[!, :MULTIPLETUMOURS] = MULTIPLETUMOUR
size(pfdt)
print(pfdt[1:5, [1,2, 3, 13, 45]])

# check how keys work in sact -> it looks like the id of patient will work ayway
print(sact_patients[1:5, :])
print(pfdt[[2569, 7736,10467, 16750, 18718],[1,27]])

# add a column for the time from diagnose to death who had surgery
TIMEDIAGNOSIStoFINAL = pfdt[:,:VITALSTATUSDATE] .- pfdt[:,:DIAGNOSISDATEBEST]
pfdt[!, :TIMEDIAGNOSITOFINALEVENT] = TIMEDIAGNOSIStoFINAL

# check the frequancy of the final outcome
freq_fs = combine(groupby(pfdt, :NEWVITALSTATUS), nrow => :Freq)
freq_fs

# exclude pateints wiht missing final outcome
delete!(pfdt, findall(ismissing.(pfdt.:NEWVITALSTATUS)))

# exclude pateints wiht lost to follow up final outcome
delete!(pfdt, findall(pfdt.:NEWVITALSTATUS.=="X"))

# check if we can add the stage
freq_stage = combine(groupby(pfdt, :STAGE_BEST), nrow => :Freq)
print(freq_stage)

# create a vector of missing values for the stage number
STAGENUM = ones(size(pfdt)[1])*99999

# for the boolean to work the missings should be character too
STAGEVEC = pfdt.:STAGE_BEST
STAGEVEC[findall(ismissing.(pfdt2.:STAGE_BEST))].="M"

# find all the STAGE 1 cancer
STAGENUM[findall(startswith.(STAGEVEC, "0"))] .= 0
STAGENUM[findall(startswith.(STAGEVEC, "1"))] .= 1
STAGENUM[findall(startswith.(STAGEVEC, "2"))] .= 2
STAGENUM[findall(startswith.(STAGEVEC, "3"))] .= 3
STAGENUM[findall(startswith.(STAGEVEC, "4"))] .= 4
STAGENUM[findall(startswith.(STAGEVEC, "5"))] .= 5
STAGENUM[findall((STAGEVEC.=="M") .|(STAGEVEC.=="6").|(STAGEVEC.=="?"
      ).|(STAGEVEC.=="U").|(STAGEVEC.=="X"))] .= -1
countmap(STAGENUM)

# add column
pfdt[!, :STAGENUMBER] = STAGENUM

# print the new summaries with new variables
pfdt_summ = describe(pfdt)
print(pfdt_summ)



# only covariates of interest
pdft_ofinterest = pfdt[:, [:PATIENTID, :SEX,:NEWVITALSTATUS, :VITALSTATUSDATE,
                           :TUMOURID, :DIAGNOSISDATEBEST, :AGE, :MULTIPLETUMOURS,
                           :TIMEDIAGNOSITOFINALEVENT, :STAGENUMBER] ]

# save the full dataset
open(string(CODEwd,"data/BDA_preprocessed.bin"), "w") do io
    serialize(io, pdft_ofinterest)
end


# y = open(deserialize, string(CODEwd,"data/BDA_preprocessed.bin"))
