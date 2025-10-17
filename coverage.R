source("compare_BCa_vs_others.R")
source("run_PEIE_addons_generic.R")
# #例：DT_Bar_WSE / B2 / COSタグ=ThiInt / DS1 / h / ut / J=2
# res <- compare_BCa_vs_others(
#   DS=1, GroupDir="DT_Bar_WSE", FileName="B2",
#   BSMethod="ThiInt", dt_code="B2",
#   thrMode="h", thrName="ut", J=2,
#   cos_tag="ThiInt",
#   alpha=0.05,
#   jk_dir="./output/DT_Bar_WSE/B2/JK-intensity",  # run_BCa_batchで保存したJK
#   make_plot=TRUE
# )    "COS","GS-","GS+","OrdSta","ThiArr","ThiInt"

# 他のCOSやJもループでOK
for (ct in c("COS","GS-","GS+","OrdSta","ThiArr","ThiInt")) {
  for (JJ in 4:5) {
    # compare_BCa_vs_others(
    #   DS=4, GroupDir="DT_Bar_WSE", FileName="B2",
    #   BSMethod=ct,               # ← ここを ct に！
    #   dt_code="B2", thrMode="h", thrName="ut", J=JJ,
    #   cos_tag=ct, make_plot=TRUE
    # )
    compare_BCa_vs_others(
      DS=1, GroupDir="DT_Ans_WSE", FileName="A3",
      BSMethod=ct,               # ← ここを ct に！
      dt_code="A3", thrMode="h", thrName="ut", J=JJ,
      cos_tag=ct, make_plot=TRUE
    )
    # compare_BCa_vs_others(
    #         DS=4, GroupDir="NDT_WSE", FileName="",
    #         BSMethod=ct,               # ← ここを ct に！
    #         dt_code="H", thrMode="h", thrName="ldt", J=JJ,
    #         cos_tag=ct, make_plot=TRUE
    # )
  }
        # for (kk in 2:6) {
        #         compare_BCa_vs_others(
        #                 DS=4, GroupDir="DT_Bar_WSE", FileName="B2",
        #                 BSMethod=ct,               # ← ここを ct に！
        #                 dt_code="B2", thrMode="h", thrName="ut", J=kk,
        #                 cos_tag=ct, make_plot=TRUE
        #         )
        #         compare_BCa_vs_others(
        #                 DS=4, GroupDir="DT_Ans_WSE", FileName="A3",
        #                 BSMethod=ct,               # ← ここを ct に！
        #                 dt_code="A3", thrMode="h", thrName="ut", J=kk,
        #                 cos_tag=ct, make_plot=TRUE
        #         )
        #         compare_BCa_vs_others(
        #                 DS=4, GroupDir="NDT_WSE", FileName="",
        #                 BSMethod=ct,               # ← ここを ct に！
        #                 dt_code="H", thrMode="h", thrName="ldt", J=kk,
        #                 cos_tag=ct, make_plot=TRUE
        #         )
        # }
}
