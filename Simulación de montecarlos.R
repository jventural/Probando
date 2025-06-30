# ----------------------------------------
# Monte Carlo EFA Simulation - BASADO EN CÓDIGO FUNCIONAL
# ----------------------------------------

# 0. Cargar librerías (igual que el código funcional)
suppressPackageStartupMessages({
  library(psych)
  library(PsyMetricTools)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(tibble)
  library(MASS)
  library(furrr)
  library(progress)
})

# 1. CONFIGURACIÓN
cores_available <- parallel::detectCores()
optimal_workers <- min(cores_available - 2, 12)
cat("🖥️  Cores disponibles:", cores_available, "| Usando:", optimal_workers, "workers\n")

plan(multisession, workers = optimal_workers)
options(future.rng.onMisuse = "ignore")
options(future.globals.maxSize = 2048^3)

# 2. FUNCIÓN EFA EXACTAMENTE COMO EL CÓDIGO FUNCIONAL

run_fa <- function(R, nfactors, method, n.obs) {
  args <- list(
    r        = R,
    nfactors = nfactors,
    rotate   = "oblimin",
    fm       = method,
    n.obs    = n.obs,
    maxit    = if (method == "ml") 10000 else NULL,  # EXACTAMENTE como el código funcional
    scores   = "tenBerge",  # puntajes más estables
    warnings = FALSE        # suprime warnings internos
  )
  # EXACTAMENTE la misma línea del código funcional
  args <- args[!vapply(args, is.null, logical(1))]
  suppressMessages(do.call(fa, args))
}

# 3. MÉTODOS EXACTAMENTE COMO EL CÓDIGO FUNCIONAL

methods <- c("ml",      # máxima verosimilitud
             "minres",  # mínimos residuos
             "pa",      # eje principal (principal axis)
             "wls",     # weighted least squares
             "gls",     # generalized least squares
             "uls")     # unweighted least squares

# 4. FUNCIÓN WLSMV USANDO PsyMetricTools (COMO EL CÓDIGO FUNCIONAL)

run_wlsmv_modern <- function(df, nfactors) {
  wlsmv_result <- tryCatch({
    
    if (ncol(df) < 3 || nrow(df) < 50 || nfactors >= ncol(df)) {
      stop("Invalid dimensions for WLSMV")
    }
    
    # Usar PsyMetricTools::EFA_modern EXACTAMENTE como el código funcional
    suppressMessages(suppressWarnings({
      PsyMetricTools::EFA_modern(
        n_factors       = nfactors,
        n_items         = ncol(df),
        name_items      = "EAI",
        data            = df,
        apply_threshold = TRUE,  # Como en el código funcional
        rotation        = "oblimin",
        estimator       = "WLSMV",
        exclude_items   = NULL
      )
    }))
    
  }, error = function(e) {
    return(structure(list(error = paste("WLSMV failed:", e$message)), class = "wlsmv_error"))
  })
  
  return(wlsmv_result)
}

# 5. SIMULACIÓN DE DATOS MEJORADA

simulate_dataset <- function(N, F, p, cat, loading_value = 0.7) {
  Lambda <- matrix(0, p, F)
  fac_assign <- rep(seq_len(F), length.out = p)
  Lambda[cbind(seq_len(p), fac_assign)] <- loading_value
  
  # Agregar variabilidad realista
  Lambda[Lambda != 0] <- Lambda[Lambda != 0] + rnorm(sum(Lambda != 0), 0, 0.05)
  Lambda[Lambda != 0] <- pmax(0.3, pmin(0.9, Lambda[Lambda != 0]))
  
  Sigma <- Lambda %*% t(Lambda) + diag(p) * (1 - loading_value^2)
  
  # Asegurar matriz definida positiva
  eigenvals <- eigen(Sigma, symmetric = TRUE)$values
  if (min(eigenvals) < 1e-6) {
    eigendecomp <- eigen(Sigma, symmetric = TRUE)
    eigenvals_corrected <- pmax(eigendecomp$values, 1e-6)
    Sigma <- eigendecomp$vectors %*% diag(eigenvals_corrected) %*% t(eigendecomp$vectors)
  }
  
  Y <- MASS::mvrnorm(N, rep(0, p), Sigma)
  cuts <- qnorm(seq(1, cat - 1) / cat)
  X_ord <- apply(Y, 2, function(x) cut(x, breaks = c(-Inf, cuts, Inf), labels = FALSE))
  
  df <- as.data.frame(X_ord)
  colnames(df) <- paste0("EAI", seq_len(p))
  list(data = df, Lambda_true = Lambda)
}

# 6. ANÁLISIS COMPLETO SIGUIENDO EL PATRÓN FUNCIONAL

analyze_dataset_functional <- function(df, Lambda_true) {
  F <- ncol(Lambda_true)
  n <- nrow(df)
  
  # 1. Matriz de correlaciones policóricas EXACTAMENTE como el código funcional
  R <- tryCatch({
    psych::polychoric(df)$rho
  }, error = function(e) {
    # Fallback: correlaciones de Spearman
    cor(df, use = "pairwise.complete.obs", method = "spearman")
  })
  
  # Asegurar que R es válida
  if (any(!is.finite(R))) {
    R[!is.finite(R)] <- 0
    diag(R) <- 1
  }
  
  # 2. Ejecutar EFA para métodos psych EXACTAMENTE como el código funcional
  efa_list <- tryCatch({
    map(methods, ~ run_fa(R = R, nfactors = F, method = .x, n.obs = n))
  }, error = function(e) {
    # Si falla, crear lista vacía
    vector("list", length(methods))
  })
  names(efa_list) <- methods
  
  # 3. Extraer métricas de métodos psych
  psych_results <- map2_dfr(efa_list, methods, function(efa, method) {
    if (is.null(efa) || inherits(efa, "try-error")) {
      return(data.frame(
        method = method, TLI = NA, CFI = NA, RMSEA = NA,
        ARB_load = NA, ARMSE_load = NA, ARB_phi = NA, ARMSE_phi = NA,
        convergence = FALSE, stringsAsFactors = FALSE
      ))
    }
    
    # Calcular métricas
    metrics <- calculate_fa_metrics_functional(efa, Lambda_true, F)
    
    data.frame(
      method = method,
      TLI = metrics$TLI, CFI = metrics$CFI, RMSEA = metrics$RMSEA,
      ARB_load = metrics$ARB_load, ARMSE_load = metrics$ARMSE_load,
      ARB_phi = metrics$ARB_phi, ARMSE_phi = metrics$ARMSE_phi,
      convergence = TRUE, stringsAsFactors = FALSE
    )
  })
  
  # 4. EFA WLSMV EXACTAMENTE como el código funcional
  wlsmv_fit <- run_wlsmv_modern(df, F)
  wlsmv_metrics <- extract_wlsmv_functional(wlsmv_fit, Lambda_true, F)
  
  # 5. Combinar resultados
  return(bind_rows(psych_results, wlsmv_metrics))
}

# 7. FUNCIONES AUXILIARES SIGUIENDO EL PATRÓN

calculate_fa_metrics_functional <- function(efa, Lambda_true, F) {
  tryCatch({
    # Cargas
    if (!is.null(efa$loadings)) {
      est_load <- as.vector(as.matrix(efa$loadings))
      true_load <- as.vector(Lambda_true)
      
      if (length(est_load) == length(true_load)) {
        valid_indices <- !is.na(est_load) & !is.na(true_load) & is.finite(est_load) & is.finite(true_load)
        if (sum(valid_indices) > 0) {
          rel_err_load <- ifelse(abs(true_load) < 1e-6, 
                                 est_load, 
                                 (est_load - true_load) / abs(true_load))
          ARB_load <- mean(rel_err_load[valid_indices], na.rm = TRUE)
          ARMSE_load <- sqrt(mean(rel_err_load[valid_indices]^2, na.rm = TRUE))
        } else {
          ARB_load <- ARMSE_load <- NA
        }
      } else {
        ARB_load <- ARMSE_load <- NA
      }
    } else {
      ARB_load <- ARMSE_load <- NA
    }
    
    # Phi (correlaciones entre factores)
    phi_true <- if (F > 1) rep(0, F*(F-1)/2) else numeric(0)
    if (F > 1 && !is.null(efa$Phi) && is.matrix(efa$Phi)) {
      est_phi <- efa$Phi[lower.tri(efa$Phi)]
      if (length(est_phi) == length(phi_true) && all(is.finite(est_phi))) {
        ARB_phi <- mean(est_phi - phi_true, na.rm = TRUE)
        ARMSE_phi <- sqrt(mean((est_phi - phi_true)^2, na.rm = TRUE))
      } else {
        ARB_phi <- ARMSE_phi <- NA
      }
    } else {
      ARB_phi <- ARMSE_phi <- NA
    }
    
    # Índices de ajuste
    TLI <- if (!is.null(efa$TLI) && is.finite(efa$TLI)) efa$TLI else NA
    CFI <- if (!is.null(efa$CFI) && is.finite(efa$CFI)) efa$CFI else NA
    RMSEA <- if (!is.null(efa$RMSEA) && length(efa$RMSEA) > 0 && is.finite(efa$RMSEA[1])) efa$RMSEA[1] else NA
    
    return(list(TLI = TLI, CFI = CFI, RMSEA = RMSEA, ARB_load = ARB_load, 
                ARMSE_load = ARMSE_load, ARB_phi = ARB_phi, ARMSE_phi = ARMSE_phi))
  }, error = function(e) {
    return(list(TLI = NA, CFI = NA, RMSEA = NA, ARB_load = NA, ARMSE_load = NA,
                ARB_phi = NA, ARMSE_phi = NA))
  })
}

# 8. EXTRACCIÓN WLSMV SIGUIENDO EL PATRÓN

extract_wlsmv_functional <- function(wlsmv_result, Lambda_true, F) {
  
  if (inherits(wlsmv_result, "wlsmv_error")) {
    return(data.frame(
      method = "wlsmv", TLI = NA, CFI = NA, RMSEA = NA,
      ARB_load = NA, ARMSE_load = NA, ARB_phi = NA, ARMSE_phi = NA,
      convergence = FALSE, stringsAsFactors = FALSE
    ))
  }
  
  tryCatch({
    # Cargas factoriales WLSMV (siguiendo el patrón del código funcional)
    wlsmv_loadings_df <- wlsmv_result$result_df %>%
      rename(item = Items) %>%
      pivot_longer(
        cols      = starts_with("f"),
        names_to  = "factor",
        values_to = "loading"
      )
    
    # Convertir a matriz para comparación
    items <- unique(wlsmv_loadings_df$item)
    factors <- paste0("f", 1:F)
    
    est_loadings_matrix <- matrix(0, length(items), F)
    for (i in 1:nrow(wlsmv_loadings_df)) {
      item_idx <- which(items == wlsmv_loadings_df$item[i])
      factor_idx <- as.numeric(gsub("f", "", wlsmv_loadings_df$factor[i]))
      if (length(item_idx) > 0 && factor_idx <= F) {
        est_loadings_matrix[item_idx, factor_idx] <- wlsmv_loadings_df$loading[i]
      }
    }
    
    # Calcular métricas de cargas
    est_load <- as.vector(est_loadings_matrix)
    true_load <- as.vector(Lambda_true)
    
    if (length(est_load) == length(true_load)) {
      valid_indices <- !is.na(est_load) & !is.na(true_load) & is.finite(est_load) & is.finite(true_load)
      if (sum(valid_indices) > 0) {
        rel_err_load <- ifelse(abs(true_load) < 1e-6, 
                               est_load, 
                               (est_load - true_load) / abs(true_load))
        ARB_load <- mean(rel_err_load[valid_indices], na.rm = TRUE)
        ARMSE_load <- sqrt(mean(rel_err_load[valid_indices]^2, na.rm = TRUE))
      } else {
        ARB_load <- ARMSE_load <- NA
      }
    } else {
      ARB_load <- ARMSE_load <- NA
    }
    
    # Índices de ajuste WLSMV (siguiendo el patrón del código funcional)
    if (!is.null(wlsmv_result$Bondades_Original)) {
      bondades <- wlsmv_result$Bondades_Original %>%
        filter(Factores == paste0("f", F))
      
      TLI <- if (nrow(bondades) > 0) as.numeric(bondades$tli.scaled[1]) else NA
      CFI <- if (nrow(bondades) > 0) as.numeric(bondades$cfi.scaled[1]) else NA
      RMSEA <- if (nrow(bondades) > 0) as.numeric(bondades$rmsea.scaled[1]) else NA
    } else {
      TLI <- CFI <- RMSEA <- NA
    }
    
    # Phi (correlaciones entre factores) - siguiendo el patrón
    if (F > 1 && !is.null(wlsmv_result$InterFactor)) {
      phi_wlsmv <- wlsmv_result$InterFactor
      est_phi <- phi_wlsmv[lower.tri(phi_wlsmv)]
      phi_true <- rep(0, F*(F-1)/2)  # Asumiendo factores ortogonales en población
      
      if (length(est_phi) == length(phi_true) && all(is.finite(est_phi))) {
        ARB_phi <- mean(est_phi - phi_true, na.rm = TRUE)
        ARMSE_phi <- sqrt(mean((est_phi - phi_true)^2, na.rm = TRUE))
      } else {
        ARB_phi <- ARMSE_phi <- NA
      }
    } else {
      ARB_phi <- ARMSE_phi <- NA
    }
    
    data.frame(
      method = "wlsmv", 
      TLI = TLI, CFI = CFI, RMSEA = RMSEA,
      ARB_load = ARB_load, ARMSE_load = ARMSE_load,
      ARB_phi = ARB_phi, ARMSE_phi = ARMSE_phi,
      convergence = TRUE, stringsAsFactors = FALSE
    )
    
  }, error = function(e) {
    data.frame(
      method = "wlsmv", TLI = NA, CFI = NA, RMSEA = NA,
      ARB_load = NA, ARMSE_load = NA, ARB_phi = NA, ARMSE_phi = NA,
      convergence = FALSE, stringsAsFactors = FALSE
    )
  })
}

# 9. CONFIGURACIÓN DE SIMULACIÓN

sample_sizes <- c(100, 200, 500)
n_factors <- 1:3
n_items <- c(5, 10)
n_cat <- c(3, 5)
n_rep <- 5  # Cambiar a 500 para análisis final

conditions <- expand_grid(N = sample_sizes, F = n_factors, p = n_items, cat = n_cat)
total_tasks <- nrow(conditions) * n_rep

all_tasks <- expand_grid(
  cond_idx = seq_len(nrow(conditions)),
  rep = seq_len(n_rep)
) %>%
  mutate(
    N = conditions$N[cond_idx],
    F = conditions$F[cond_idx], 
    p = conditions$p[cond_idx],
    cat = conditions$cat[cond_idx]
  ) %>%
  dplyr::select(-cond_idx)

chunk_size <- max(2, floor(total_tasks / (optimal_workers * 15)))
chunks <- split(all_tasks, ceiling(seq_len(nrow(all_tasks)) / chunk_size))

cat("🖥️ Total tareas:", total_tasks, "| Chunks:", length(chunks), "\n")
cat("⚡ Métodos: ML, MINRES, PA, WLS, GLS, ULS + WLSMV\n")
cat("🎯 Basado en código funcional probado\n")

# 10. PROCESAMIENTO SIGUIENDO EL PATRÓN FUNCIONAL

process_chunk_functional <- function(chunk_df) {
  chunk_results <- vector("list", nrow(chunk_df))
  
  for (i in seq_len(nrow(chunk_df))) {
    task <- chunk_df[i, ]
    
    result <- tryCatch({
      set.seed(task$rep * 1000 + task$N + task$F * 10 + task$p)
      sim <- simulate_dataset(task$N, task$F, task$p, task$cat)
      analyze_dataset_functional(sim$data, sim$Lambda_true)
    }, error = function(e) {
      # Resultado de emergencia
      methods_all <- c("ml", "minres", "pa", "wls", "gls", "uls", "wlsmv")
      data.frame(
        method = methods_all, TLI = NA, CFI = NA, RMSEA = NA, 
        ARB_load = NA, ARMSE_load = NA, ARB_phi = NA, ARMSE_phi = NA, 
        convergence = FALSE, stringsAsFactors = FALSE
      )
    })
    
    result$N <- task$N
    result$F <- task$F
    result$p <- task$p
    result$cat <- task$cat
    result$rep <- task$rep
    chunk_results[[i]] <- result
  }
  
  return(do.call(rbind, chunk_results))
}

# 11. EJECUCIÓN

pb <- progress_bar$new(
  format = "[:bar] :percent | :current/:total | ⏱️ :elapsed | ETA: :eta",
  total = length(chunks), clear = FALSE, width = 80
)

cat("\n🚀 SIMULACIÓN BASADA EN CÓDIGO FUNCIONAL\n")
cat("📊 Siguiendo exactamente el patrón que ya funciona\n")
cat("⚡ run_fa() + PsyMetricTools::EFA_modern()\n")
cat("🛡️ Mismos parámetros y estructura\n\n")

start_time <- Sys.time()
chunk_results <- vector("list", length(chunks))

for (i in seq_along(chunks)) {
  cat("📊 Chunk", i, "/", length(chunks), "(", nrow(chunks[[i]]), "tareas)...")
  
  chunk_start <- Sys.time()
  chunk_results[[i]] <- process_chunk_functional(chunks[[i]])
  chunk_time <- as.numeric(Sys.time() - chunk_start, units = "secs")
  
  pb$tick()
  cat(" ✅", round(chunk_time, 1), "s\n")
  
  if (i %% max(1, length(chunks) %/% 4) == 0 || i == length(chunks)) {
    elapsed <- as.numeric(Sys.time() - start_time, units = "secs")
    rate <- i / elapsed * 60
    cat("🚀 Velocidad:", round(rate, 1), "chunks/min\n\n")
  }
  
  flush.console()
}

# 12. RESULTADOS FINALES

results_sim <- do.call(rbind, chunk_results)
end_time <- Sys.time()
total_time <- end_time - start_time

cat("\n🎉 ¡SIMULACIÓN COMPLETADA!\n")
cat("⏱️ Tiempo:", round(total_time, 2), attr(total_time, "units"), "\n")

# Estadísticas por método
method_stats <- results_sim %>%
  group_by(method) %>%
  summarise(
    total_cases = n(),
    convergence_rate = mean(convergence, na.rm = TRUE) * 100,
    success_rate = mean(!is.na(TLI), na.rm = TRUE) * 100,
    .groups = "drop"
  )

cat("\n📊 ESTADÍSTICAS POR MÉTODO:\n")
print(method_stats)

# Resumen principal
summary_sim <- results_sim %>%
  filter(!is.na(TLI)) %>%
  group_by(N, F, p, cat, method) %>%
  summarise(
    n_valid = n(),
    TLI_mean = mean(TLI, na.rm = TRUE),
    CFI_mean = mean(CFI, na.rm = TRUE), 
    RMSEA_mean = mean(RMSEA, na.rm = TRUE),
    ARB_load_mean = mean(ARB_load, na.rm = TRUE),
    ARMSE_load_mean = mean(ARMSE_load, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n📈 RESUMEN DE RESULTADOS:\n")
print(head(summary_sim, 20))

plan(sequential)
cat("\n🏁 Simulación basada en código funcional completada.\n")
cat("✅ Misma estructura que el análisis que ya funciona\n")

# Proyección para 500 réplicas
cat("\n🔮 PARA 500 RÉPLICAS:\n")
tiempo_estimado_500 <- as.numeric(total_time) * (500 / n_rep)
cat("⏰ Tiempo estimado:", round(tiempo_estimado_500 / 60, 1), "minutos\n")
cat("✅ Código funcional escalado a Monte Carlo\n")