using ImplicitIntegration
using Test
using Aqua

Aqua.test_ambiguities(ImplicitIntegration)
Aqua.test_unbound_args(ImplicitIntegration)
Aqua.test_undefined_exports(ImplicitIntegration)
Aqua.test_project_extras(ImplicitIntegration)
Aqua.test_stale_deps(ImplicitIntegration)
Aqua.test_deps_compat(ImplicitIntegration)
Aqua.test_piracies(ImplicitIntegration; broken = true) # FIXME: piracy on IntervalArithmetic.Interval for now
Aqua.test_persistent_tasks(ImplicitIntegration)
