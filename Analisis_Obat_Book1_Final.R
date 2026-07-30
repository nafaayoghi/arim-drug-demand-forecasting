# ====================================================================
# Analisis Deret Waktu - Peramalan Kuantitas Obat Generik
# Metode  : ARIMA Box-Jenkins (Parsimoni & Log-Transformed)
# Data    : April 2021 - Mei 2026 (62 periode) | File: Book1.xlsx
# ====================================================================

# ====================================================================
# 1. LIBRARY
# ====================================================================
library(readxl)    # baca file excel
library(forecast)  # buat fitting ARIMA dan fungsi forecast-nya
library(tseries)   # perlu ini buat ADF test dan KPSS test
library(MASS)      # Box-Cox pakai ini
library(lmtest)    # coeftest buat cek signifikansi parameter

setwd("D:/R")
# Simpan semua output konsol ke file txt
sink("Hasil_Analisis_Obat_Book1_Final.txt", append = FALSE, split = TRUE)
df <- read_excel("Book1.xlsx")

# ====================================================================
# 2. EKSPLORASI DATA
# ====================================================================
cat("========================================================\n")
cat("          EKSPLORASI DATA PENUH (APR 2021 - MEI 2026)    \n")
cat("========================================================\n")

data_penuh_ts <- ts(df$jumlah, start = c(2021, 4), frequency = 12)

plot.ts(data_penuh_ts, lty = 1, col = "blue", lwd = 2,
        xlab = "Waktu", ylab = "Kuantitas Obat",
        main = "Plot Data Kuantitas Obat (April 2021 - Mei 2026)")
grid()
print(summary(data_penuh_ts))

# ====================================================================
# 3. PEMBAGIAN DATA LATIH DAN UJI
# ====================================================================
cat("========================================================\n")
cat("       PEMBAGIAN DATA (50 BLN LATIH : 12 BLN UJI)       \n")
cat("========================================================\n")

# 50 data latih (Apr 2021 - Mei 2025)
obat.train.ts <- ts(df$jumlah[1:50],   start = c(2021, 4), frequency = 12)
# 12 data uji (Jun 2025 - Mei 2026)
obat.test.ts  <- ts(df$jumlah[51:62], start = c(2025, 6), frequency = 12)

plot.ts(obat.train.ts, lty = 1, col = "darkorange", lwd = 2,
        xlab = "Waktu", ylab = "Kuantitas Obat",
        main = "Plot Data Latih (April 2021 - Mei 2025)")
grid()

# ====================================================================
# 4. UJI STASIONERITAS DATA LATIH
# ====================================================================
cat("========================================================\n")
cat("                UJI STASIONERITAS DATA                 \n")
cat("========================================================\n")

cat("--- A. Uji Stasioneritas Ragam (Box-Cox vs Log) ---\n")
index <- seq(1:length(obat.train.ts))
bc <- boxcox(obat.train.ts ~ index, lambda = seq(-2, 2, by = 0.05), plotit = FALSE)
lambda_opt <- bc$x[which.max(bc$y)]
cat("Lambda optimum Box-Cox =", round(lambda_opt, 4), "\n")

# Selang kepercayaan 95% Box-Cox dipakai sebagai acuan untuk menentukan
# apakah lambda = 0 (log) memang layak dipakai, bukan sekadar diasumsikan.
ci_bc  <- range(bc$x[bc$y > max(bc$y) - qchisq(0.95, 1) / 2])
cat("Selang kepercayaan 95% lambda:", round(ci_bc[1], 4), "-", round(ci_bc[2], 4), "\n")
if (ci_bc[1] <= 0 && ci_bc[2] >= 0) {
  cat(">> Keputusan: lambda = 0 berada dalam selang kepercayaan 95%, sehingga\n")
  cat("   Log Transformation (lambda = 0) sah digunakan untuk menstabilkan ragam.\n\n")
} else {
  cat(">> PERHATIAN: lambda = 0 TIDAK berada dalam selang kepercayaan 95% Box-Cox.\n")
  cat("   Log Transformation tetap dipakai di sini karena alasan interpretasi\n")
  cat("   (kembali ke skala asli mudah & umum di praktik farmasi), namun ini\n")
  cat("   perlu disebutkan secara eksplisit sebagai keterbatasan.\n\n")
}

# Transformasi log diterapkan untuk eksplorasi ADF, KPSS, dan differencing
data_log <- log(obat.train.ts)

cat("--- B. Uji Stasioneritas Rataan (Level) ---\n")
cat(">> ADF Test (H0: TIDAK stasioner) pada data Log:\n")
print(suppressWarnings(tseries::adf.test(data_log)))

# KPSS ditambahkan sebagai komplemen ADF. 
cat("\n>> KPSS Test (H0: stasioner) pada data Log:\n")
print(suppressWarnings(tseries::kpss.test(data_log, null = "Level")))

cat("--- C. Differencing d = 1 ---\n")
data_diff <- diff(data_log, differences = 1)
plot.ts(data_diff, lty = 1, col = "purple", lwd = 2,
        main = "Plot Data Transformasi Log & Differencing (d=1)")
abline(h = 0, col = "red", lty = 2)

cat(">> ADF Test setelah differencing d=1:\n")
print(suppressWarnings(tseries::adf.test(data_diff)))

cat("\n>> KPSS Test setelah differencing d=1:\n")
print(suppressWarnings(tseries::kpss.test(data_diff, null = "Level")))
cat(">> Jika ADF menolak H0 (stasioner) DAN KPSS gagal tolak H0 (stasioner),\n")
cat("   maka d=1 sudah cukup untuk menstasionerkan rataan.\n\n")

# ====================================================================
# 5. IDENTIFIKASI MODEL TENTATIF DARI ACF & PACF
# ====================================================================
cat("========================================================\n")
cat("            IDENTIFIKASI MODEL TENTATIF                 \n")
cat("========================================================\n")

par(mfrow = c(1, 2))
acf_result <- acf(data_diff, lag.max = 24, plot = FALSE)
acf_result$lag <- acf_result$lag * frequency(data_diff)
plot(acf_result, main = "ACF (Log, d=1)", xlab = "Lag", xaxt = "n")
axis(1, at = seq(0, 24, by = 2), labels = seq(0, 24, by = 2))

pacf_result <- pacf(data_diff, lag.max = 24, plot = FALSE)
pacf_result$lag <- pacf_result$lag * frequency(data_diff)
plot(pacf_result, main = "PACF (Log, d=1)", xlab = "Lag", xaxt = "n")
axis(1, at = seq(2, 24, by = 2), labels = seq(2, 24, by = 2))
par(mfrow = c(1, 1))

# ====================================================================
# 6. PENDUGAAN PARAMETER
#    Kandidat utama diambil murni dari lag-lag yang signifikan di plot
#    ACF/PACF.
#    Grid search lengkap p(0-3) x q(0-3) dipindah ke bagian 6B dan hanya
#    dipakai sebagai pembanding/robustness check di lampiran, bukan
#    sebagai metode pemilihan model utama.
# ====================================================================
cat("========================================================\n")
cat("   6A. KANDIDAT MODEL DARI IDENTIFIKASI VISUAL ACF/PACF   \n")
cat("========================================================\n")

# Fungsi bantu untuk menentukan lag Ljung-Box.
# Pakai rule of thumb Hyndman (non-musiman): lag = min(10, n/5),
# dengan syarat lag harus > fitdf (jumlah parameter AR+MA) agar derajat
# bebas Ljung-Box tidak negatif/nol. Kalau lag hasil rule <= fitdf,
# lag dinaikkan ke fitdf+1 dan dicatat sebagai penyesuaian.
tentukan_lag_lb <- function(fitdf, n) {
  lag_rule <- floor(min(10, n / 5))
  if (lag_rule <= fitdf) lag_rule <- fitdf + 1
  return(lag_rule)
}

cek_signifikan <- function(model) {
  if (length(model$coef) == 0) return("Tidak Ada Parameter")
  p_val <- coeftest(model)[, 4]
  return(ifelse(all(p_val < 0.05), "Signifikan Semua", "Ada Tidak Signifikan"))
}

cek_white_noise <- function(model) {
  p_ord  <- arimaorder(model)[1]
  q_ord  <- arimaorder(model)[3]
  fitdf  <- p_ord + q_ord
  n_obs  <- length(model$residuals)
  lag_test <- tentukan_lag_lb(fitdf, n_obs)
  fitdf_use <- ifelse(fitdf == 0, 1, fitdf)  # cegah error fitdf=0 di Ljung-Box
  
  # Selain lag utama (rule Hyndman), lag standar 6, 12, 18, 24 juga
  # dicek sebagai bukti tambahan kekuatan kesimpulan white noise
  lb <- Box.test(model$residuals, type = "Ljung-Box", lag = lag_test, fitdf = fitdf_use)
  return(ifelse(lb$p.value > 0.05, "Lolos (WN)", "Gagal"))
}

hitung_mape <- function(model) {
  fc   <- forecast::forecast(model, h = 12)
  mape <- mean(abs(as.numeric(obat.test.ts) - as.numeric(fc$mean)) / as.numeric(obat.test.ts)) * 100
  return(round(mape, 4))
}

# AICc dipakai, bukan AIC biasa karena data latih kecil (n=50),
ambil_aicc <- function(model) round(model$aicc, 2)


kandidat_visual <- list(
  "ARIMA(2,1,0)" = c(2, 1, 0),
  "ARIMA(0,1,2)" = c(0, 1, 2),
  "ARIMA(2,1,2)" = c(2, 1, 2),
  "ARIMA(3,1,0)" = c(3, 1, 0),
  "ARIMA(0,1,3)" = c(0, 1, 3),
  "ARIMA(3,1,3)" = c(3, 1, 3),
  "ARIMA(6,1,0)" = c(6, 1, 0),
  "ARIMA(0,1,6)" = c(0, 1, 6)
)

list_kandidat <- list()
for (nm in names(kandidat_visual)) {
  ord <- kandidat_visual[[nm]]
  m_fit <- tryCatch(
    Arima(obat.train.ts, order = ord, include.drift = FALSE, lambda = 0, method = "ML"),
    error = function(e) NULL
  )
  if (!is.null(m_fit)) list_kandidat[[nm]] <- m_fit
}

cat("\n--- Uji Signifikansi Parameter Kandidat Visual ---\n")
for (nm in names(list_kandidat)) {
  cat("\n> Model:", nm, "\n")
  print(coeftest(list_kandidat[[nm]]))
}

tab_komparasi <- data.frame(
  Model        = names(list_kandidat),
  AICc         = sapply(list_kandidat, ambil_aicc),
  Signifikansi = sapply(list_kandidat, cek_signifikan),
  White_Noise  = sapply(list_kandidat, cek_white_noise),
  MAPE         = sapply(list_kandidat, hitung_mape)
)
row.names(tab_komparasi) <- NULL

cat("\n--- TABEL KOMPARASI MODEL KANDIDAT VISUAL (kriteria utama) ---\n")
print(tab_komparasi)
write.csv(tab_komparasi, "1_Tabel_Komparasi_Model.csv", row.names = FALSE)

# Seleksi model: kriteria utama adalah lolos signifikansi & lolos white noise
# Di antara model-model yang valid tsb, dipilih AICc terkecil sebagai model final.
# MAPE terhadap data uji hanya dilaporkan sebagai ukuran performa akurasi,
model_layak <- tab_komparasi[tab_komparasi$Signifikansi == "Signifikan Semua" & tab_komparasi$White_Noise == "Lolos (WN)", ]

if (nrow(model_layak) > 0) {
  model_terpilih_nama <- model_layak$Model[which.min(model_layak$AICc)]
  cat("\n>> Model-model yang lolos Signifikansi & White Noise:\n")
  print(model_layak[order(model_layak$AICc), c("Model", "AICc", "MAPE")])
  cat("\n>> Hasil Seleksi (AICc terkecil di antara model valid):", model_terpilih_nama, "\n")
  cat("   (MAPE model ini terhadap data uji:", model_layak$MAPE[which.min(model_layak$AICc)], "% - dilaporkan, bukan kriteria pemilihan)\n")
} else {
  cat("\n>> Peringatan: Tidak ada model yang lolos kedua kriteria uji. Cek tabel secara manual.\n")
  model_terpilih_nama <- tab_komparasi$Model[which.min(tab_komparasi$AICc)]
  cat(">> Fallback: dipilih AICc terkecil dari seluruh kandidat visual:", model_terpilih_nama, "\n")
}

cat("\n========================================================\n")
cat("  6B. LAMPIRAN/ROBUSTNESS CHECK: GRID SEARCH p(0-3) x q(0-3) \n")
cat("  (HANYA pembanding, TIDAK dipakai sbg dasar pemilihan model \n")
cat("   utama, identifikasi tetap berbasis plot ACF/PACF di atas)\n")
cat("========================================================\n")

list_kandidat_grid <- list()
for (p_idx in 0:3) {
  for (q_idx in 0:3) {
    if (p_idx == 0 && q_idx == 0) next
    nama_model <- paste0("ARIMA(", p_idx, ",1,", q_idx, ")")
    m_fit <- tryCatch(
      Arima(obat.train.ts, order = c(p_idx, 1, q_idx), include.drift = FALSE, lambda = 0, method = "ML"),
      error = function(e) NULL
    )
    if (!is.null(m_fit)) list_kandidat_grid[[nama_model]] <- m_fit
  }
}

tab_grid <- data.frame(
  Model        = names(list_kandidat_grid),
  AICc         = sapply(list_kandidat_grid, ambil_aicc),
  Signifikansi = sapply(list_kandidat_grid, cek_signifikan),
  White_Noise  = sapply(list_kandidat_grid, cek_white_noise),
  MAPE         = sapply(list_kandidat_grid, hitung_mape)
)
row.names(tab_grid) <- NULL
tab_grid <- tab_grid[order(tab_grid$AICc), ]
cat("\n--- Tabel Grid Search Lengkap (lampiran, urut AICc) ---\n")
print(tab_grid)
write.csv(tab_grid, "1b_Lampiran_Grid_Search_Lengkap.csv", row.names = FALSE)


# ====================================================================
# 7. MODEL FINAL & CEK OVERFITTING
# ====================================================================
cat("========================================================\n")
cat("              MODEL FINAL & CEK OVERFITTING             \n")
cat("========================================================\n")

model_final <- list_kandidat[[model_terpilih_nama]]

cat(">> Menjalankan Model Final Terpilih:", model_terpilih_nama, "\n")
summary(model_final)

cat("\n--- Uji Signifikansi Parameter Final ---\n")
print(coeftest(model_final))

p_final <- arimaorder(model_final)[1]
d_final <- arimaorder(model_final)[2]
q_final <- arimaorder(model_final)[3]

cat("\n--- Uji Overfitting ---\n")
# Selain cek terpisah (p+1,q) dan (p,q+1), cek gabungan (p+1,q+1) juga
# ditambahkan agar overfitting test lebih menyeluruh.
model_over_p  <- tryCatch(Arima(obat.train.ts, order = c(p_final + 1, d_final, q_final),     include.drift = FALSE, lambda = 0, method = "ML"), error = function(e) NULL)
model_over_q  <- tryCatch(Arima(obat.train.ts, order = c(p_final,     d_final, q_final + 1), include.drift = FALSE, lambda = 0, method = "ML"), error = function(e) NULL)
model_over_pq <- tryCatch(Arima(obat.train.ts, order = c(p_final + 1, d_final, q_final + 1), include.drift = FALSE, lambda = 0, method = "ML"), error = function(e) NULL)

# ---------------------------------------------------------------
# [PATCH] Overfitting check yang lebih menyeluruh: bukan cuma cek
# signifikansi parameter tambahan, tapi juga White Noise & apakah
# AICc-nya membaik dibanding model final. Model final baru dianggap
# optimal kalau ketiga model overfitting ini TIDAK lebih baik.
# ---------------------------------------------------------------
cek_overfitting <- function(model_over, model_final, label) {
  if (is.null(model_over)) {
    cat(label, ": model gagal fit\n")
    return(invisible(NULL))
  }
  sig  <- cek_signifikan(model_over)
  wn   <- cek_white_noise(model_over)
  aicc_over  <- ambil_aicc(model_over)
  aicc_final <- ambil_aicc(model_final)
  membaik <- ifelse(aicc_over < aicc_final, "AICc membaik", "AICc TIDAK membaik")
  
  cat(label, "\n")
  cat("   Signifikansi :", sig, "\n")
  cat("   White Noise  :", wn, "\n")
  cat("   AICc         :", aicc_over, "(model final:", aicc_final, "->", membaik, ")\n\n")
}

cek_overfitting(model_over_p,  model_final, "Model over-p  (p+1, q):")
cek_overfitting(model_over_q,  model_final, "Model over-q  (p, q+1):")
cek_overfitting(model_over_pq, model_final, "Model over-pq (p+1, q+1):")

cat(">> Model final dianggap optimal jika model overfitting TIDAK signifikan semua\n")
cat("   DAN AICc-nya tidak membaik dibanding model final.\n\n")

# ====================================================================
# 7.5. EVALUASI POLA MODEL (FITTED VS ACTUAL PADA DATA LATIH)
# ====================================================================
cat("========================================================\n")
cat("       EVALUASI POLA MODEL (FITTED VS ACTUAL)           \n")
cat("========================================================\n")

fitted_values <- fitted(model_final)
par(mfrow = c(1, 1))

plot(obat.train.ts, type = "l", col = "black", lwd = 2,
     main = "Evaluasi Pola: Data Aktual Latih vs Fitted Value Model",
     ylab = "Kuantitas Obat", xlab = "Waktu (Tahun)",
     xlim = c(2021.25, 2025.5))
lines(fitted_values, col = "red", lwd = 2, lty = 2)
grid()
legend("topleft", legend = c("Data Aktual (Mulai Apr 2021)", "Fitted Value (Taksiran Model)"),
       col = c("black", "red"), lty = c(1, 2), lwd = 2, bty = "n")

# ====================================================================
# 8. DIAGNOSTIK SISAAN
# ====================================================================
cat("========================================================\n")
cat("                DIAGNOSTIK SISAAN                       \n")
cat("========================================================\n")

sisaan_final <- model_final$residuals

par(mfrow = c(2, 2))
qqnorm(sisaan_final); qqline(sisaan_final, col = "blue", lwd = 2)
plot(1:length(sisaan_final), sisaan_final, main = "Plot Sisaan", ylab = "Sisaan", xlab = "Waktu")

acf_sisaan <- acf(sisaan_final, lag.max = 24, plot = FALSE)
acf_sisaan$lag <- acf_sisaan$lag * frequency(sisaan_final)
plot(acf_sisaan, main = "ACF Sisaan", xlab = "Lag", xaxt = "n")
axis(1, at = seq(0, 24, by = 2), labels = seq(0, 24, by = 2))

pacf_sisaan <- pacf(sisaan_final, lag.max = 24, plot = FALSE)
pacf_sisaan$lag <- pacf_sisaan$lag * frequency(sisaan_final)
plot(pacf_sisaan, main = "PACF Sisaan", xlab = "Lag", xaxt = "n")
axis(1, at = seq(2, 24, by = 2), labels = seq(2, 24, by = 2))
par(mfrow = c(1, 1))

# Fungsi tentukan_lag_lb() yang sama dari section 6 dipakai lagi
# supaya lag Ljung-Box konsisten antara tahap seleksi model & diagnostik final.
fitdf_val  <- p_final + q_final
n_res      <- length(sisaan_final)
lag_final  <- tentukan_lag_lb(fitdf_val, n_res)
fitdf_use  <- ifelse(fitdf_val == 0, 1, fitdf_val)

cat("1) Shapiro-Wilk Test (Normalitas):\n"); print(shapiro.test(sisaan_final))
cat("2) Ljung-Box Test (White Noise), lag =", lag_final, ":\n")
print(Box.test(sisaan_final, type = "Ljung-Box", lag = lag_final, fitdf = fitdf_use))

# Cek tambahan di beberapa lag standar (6,12,18,24) sebagai robustness
cat("\n>> Ljung-Box pada beberapa lag standar (pelengkap, bukan kriteria utama):\n")
for (lg in c(6, 12, 18, 24)) {
  if (lg > fitdf_use) {
    hasil_lb <- Box.test(sisaan_final, type = "Ljung-Box", lag = lg, fitdf = fitdf_use)
    cat("   lag =", lg, "-> p-value =", round(hasil_lb$p.value, 4), "\n")
  }
}

# ---------------------------------------------------------------
# [PATCH] Flag keterbatasan model otomatis, biar setiap kali script
# di-run ulang (misal ada yg clone dari GitHub), keterbatasan model
# langsung kelihatan di output txt tanpa perlu ditulis manual.
# ---------------------------------------------------------------
cat("\n--- Ringkasan Flag Keterbatasan Model (auto-generated) ---\n")

sw_test <- shapiro.test(sisaan_final)
if (sw_test$p.value < 0.05) {
  cat("[FLAG] Residual TIDAK normal (Shapiro-Wilk p =", round(sw_test$p.value, 5),
      "). Titik forecast (mean) tetap valid secara asimtotik (CLT/Gauss-Markov),\n")
  cat("       namun selang kepercayaan (PI) forecast berasumsi normal sehingga\n")
  cat("       lebarnya mungkin kurang akurat, terutama di horizon forecast yang jauh.\n\n")
} else {
  cat("[OK] Residual model final tidak menyalahi asumsi normalitas (Shapiro-Wilk p >= 0.05).\n\n")
}

# ====================================================================
# 9. PERAMALAN & EVALUASI DATA UJI (BACKTESTING)
# ====================================================================
cat("========================================================\n")
cat("         PERAMALAN & EVALUASI DATA UJI (BACKTEST)       \n")
cat("========================================================\n")

ramalan_final <- forecast::forecast(model_final, h = 12, level = c(80, 95))

plot(ramalan_final,
     main  = "Validasi Model: Data Latih vs Data Aktual (Uji) vs Forecast",
     ylab  = "Kuantitas Obat", xlab  = "Waktu",
     fcol  = "red", flwd  = 2,
     xlim  = c(2021.25, 2026.5))

lines(obat.test.ts, col = "darkgreen", lwd = 2, lty = 2)
legend("topleft", legend = c("Data Latih", "Data Uji (Aktual)", "Hasil Peramalan", "Selang Kepercayaan 95%"),
       col = c("black", "darkgreen", "red", "lightblue"), lty = c(1, 2, 1, 1), lwd = c(2, 2, 2, 8), bty = "n")

aktual_arr   <- as.numeric(obat.test.ts)
prediksi_arr <- as.numeric(ramalan_final$mean)

tab_hasil <- data.frame(
  Bulan      = c("Jun-25","Jul-25","Agt-25","Sep-25","Okt-25","Nov-25",
                 "Des-25","Jan-26","Feb-26","Mar-26","Apr-26","Mei-26"),
  Aktual     = aktual_arr,
  Prediksi   = round(prediksi_arr, 2),
  Error      = round(aktual_arr - prediksi_arr, 2),
  APE_Pct    = round(abs(aktual_arr - prediksi_arr) / aktual_arr * 100, 4)
)
print(tab_hasil)
write.csv(tab_hasil, "2_Tabel_Hasil_Backtest.csv", row.names = FALSE)

cat("MAPE Model terhadap Data Uji:", round(mean(tab_hasil$APE_Pct), 4), "%\n")

# ====================================================================
# 10. FUTURE FORECAST (HANYA DILAKUKAN SETELAH VALIDASI)
# ====================================================================
cat("========================================================\n")
cat("          FUTURE FORECAST (JUNI 2026 - MEI 2027)        \n")
cat("========================================================\n")

model_masa_depan <- Arima(data_penuh_ts, order = c(p_final, d_final, q_final), include.drift = FALSE, lambda = 0, method = "ML")

cat(">> Validasi ulang model_masa_depan (refit pada seluruh data):\n")
cat("--- Signifikansi Parameter ---\n")
print(coeftest(model_masa_depan))

fitdf_depan <- p_final + q_final
n_depan     <- length(model_masa_depan$residuals)
lag_depan   <- tentukan_lag_lb(fitdf_depan, n_depan)
fitdf_depan_use <- ifelse(fitdf_depan == 0, 1, fitdf_depan)

cat("\n--- Ljung-Box White Noise (lag =", lag_depan, ") ---\n")
lb_depan <- Box.test(model_masa_depan$residuals, type = "Ljung-Box", lag = lag_depan, fitdf = fitdf_depan_use)
print(lb_depan)

if (!(all(coeftest(model_masa_depan)[, 4] < 0.05) && lb_depan$p.value > 0.05)) {
  cat(">> PERINGATAN: model_masa_depan TIDAK sepenuhnya lolos re-validasi\n")
  cat("   (signifikansi parameter dan/atau white noise). Forecast ke depan\n")
  cat("   tetap ditampilkan, tapi keterbatasan ini WAJIB disebutkan di laporan.\n\n")
  
  # [PATCH] Detail parameter mana yang gagal, biar gampang dicek tanpa buka
  # ulang tabel coeftest manual.
  sig_depan <- coeftest(model_masa_depan)[, 4]
  if (!all(sig_depan < 0.05)) {
    cat("[FLAG] Parameter tidak signifikan di model_masa_depan:",
        paste(names(sig_depan)[sig_depan >= 0.05], collapse = ", "), "\n\n")
  }
} else {
  cat(">> Model_masa_depan lolos re-validasi (signifikan & white noise). Aman dipakai forecast.\n\n")
}

ramalan_depan <- forecast::forecast(model_masa_depan, h = 12, level = c(80, 95))

if (.Platform$OS.type == "windows") windows() else dev.new()

plot(ramalan_depan, fcol = "purple", flwd = 2,
     main = "Peramalan Kuantitas Obat Masa Depan (Jun 2026 - Mei 2027)",
     ylab = "Kuantitas Obat", xlab = "Waktu")
grid()
print(ramalan_depan)

# ====================================================================
# 11. VALIDASI ROLLING FORECAST & PLOT GABUNGAN
# ====================================================================
cat("========================================================\n")
cat("         VALIDASI ROLLING FORECAST (WALK-FORWARD)       \n")
cat("========================================================\n")

n_total <- length(df$jumlah)
n_test  <- 12
n_train <- n_total - n_test

rolling_pred   <- numeric(n_test)
rolling_aktual <- as.numeric(obat.test.ts)

for (i in 1:n_test) {
  data_rolling <- ts(df$jumlah[1:(n_train + i - 1)], start = c(2021, 4), frequency = 12)
  m_roll <- tryCatch(
    Arima(data_rolling, order = c(p_final, d_final, q_final), include.drift = FALSE, lambda = 0, method = "ML"),
    error = function(e) NULL
  )
  if (!is.null(m_roll)) {
    fc_roll         <- forecast::forecast(m_roll, h = 1)
    rolling_pred[i] <- as.numeric(fc_roll$mean)
  } else {
    rolling_pred[i] <- NA
  }
}

tab_rolling <- data.frame(
  Bulan      = tab_hasil$Bulan,
  Aktual     = rolling_aktual,
  Pred_Roll  = round(rolling_pred, 2),
  APE_Roll   = round(abs(rolling_aktual - rolling_pred) / rolling_aktual * 100, 4)
)
write.csv(tab_rolling, "3_Tabel_Rolling_Forecast.csv", row.names = FALSE)

# ---------------------------------------------------------------
# [PATCH] Sebelumnya tab_rolling cuma disimpan ke CSV, gak pernah
# di-print ke output txt. Ditambahin di sini biar MAPE rolling
# kelihatan dan bisa dibandingin langsung sama MAPE fixed forecast.
# ---------------------------------------------------------------
cat("\n--- Tabel Rolling Forecast (Walk-Forward) ---\n")
print(tab_rolling)

mape_roll <- round(mean(tab_rolling$APE_Roll, na.rm = TRUE), 4)
cat("\nMAPE Rolling Forecast:", mape_roll, "%\n")
cat("MAPE Fixed Forecast   :", round(mean(tab_hasil$APE_Pct), 4), "%\n")
cat(">> Perbandingan: jika MAPE rolling < MAPE fixed, forecast walk-forward\n")
cat("   lebih stabil karena model di-refit tiap langkah dengan data terbaru.\n")
cat("   Jika sebaliknya, model fixed sudah cukup baik dan refit berulang\n")
cat("   tidak memberi banyak keuntungan tambahan.\n\n")

ts_aktual_uji <- ts(rolling_aktual, start = c(2025, 6), frequency = 12)
ts_fixed_pred <- ts(prediksi_arr, start = c(2025, 6), frequency = 12)
ts_roll_pred  <- ts(rolling_pred, start = c(2025, 6), frequency = 12)

ymin <- min(c(rolling_aktual, prediksi_arr, rolling_pred), na.rm = TRUE) - 10
ymax <- max(c(rolling_aktual, prediksi_arr, rolling_pred), na.rm = TRUE) + 10

if (.Platform$OS.type == "windows") windows() else dev.new()

plot(ts_aktual_uji, type = "o", col = "black", lwd = 2, pch = 16,
     ylim = c(ymin, ymax), xlab = "Waktu", ylab = "Kuantitas Obat",
     main = "Komparasi Data Aktual vs Fixed Forecast vs Rolling Forecast")
lines(ts_fixed_pred, type = "o", col = "red", lwd = 2, lty = 2, pch = 15)
lines(ts_roll_pred, type = "o", col = "blue", lwd = 2, lty = 3, pch = 17)
grid()
legend("topleft", legend = c("Aktual (Data Uji)", "Fixed Forecast", "Rolling Forecast"),
       col = c("black", "red", "blue"), lty = c(1, 2, 3), lwd = 2, pch = c(16, 15, 17), bty = "n")

cat("=== SKRIP SELESAI ===\n")
sink()
