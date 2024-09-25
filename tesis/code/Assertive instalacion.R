# Lista de dependencias de assertive en Bitbucket
dependencies <- c(
  "assertive.types", "assertive.numbers", "assertive.strings", "assertive.datetimes",
  "assertive.files", "assertive.sets", "assertive.matrices", "assertive.models",
  "assertive.data", "assertive.data.uk", "assertive.data.us", "assertive.reflection",
  "assertive.code"
)

# Función para instalar paquetes desde Bitbucket
install_bitbucket_dependency <- function(package_name) {
  message(paste("Instalando", package_name, "desde Bitbucket"))
  tryCatch({
    devtools::install_bitbucket(paste0("richierocks/", package_name))
  }, warning = function(w) {
    message(paste("Advertencia durante la instalación de", package_name, ":", w))
  }, error = function(e) {
    message(paste("Error al instalar", package_name, ":", e))
  })
}

# Intentar instalar cada dependencia
for (dep in dependencies) {
  install_bitbucket_dependency(dep)
}

# Intentar instalar assertive después de las dependencias
message("Intentando instalar el paquete 'assertive' desde Bitbucket")
tryCatch({
  devtools::install_bitbucket("richierocks/assertive")
}, warning = function(w) {
  message("Advertencia durante la instalación de 'assertive':", w)
}, error = function(e) {
  message("Error durante la instalación de 'assertive':", e)
})
