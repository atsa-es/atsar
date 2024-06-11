// Generated by rstantools.  Do not edit by hand.

/*
    atsarpackage is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    atsarpackage is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with atsarpackage.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#ifndef USE_STANC3
#define USE_STANC3
#endif
#include <rstan/rstaninc.hpp>
// Code generated by stanc v2.32.2
#include <stan/model/model_header.hpp>
namespace model_rw_namespace {
using stan::model::model_base_crtp;
using namespace stan::math;
stan::math::profile_map profiles__;
static constexpr std::array<const char*, 38> locations_array__ =
  {" (found before start of program)",
  " (in 'rw', line 8, column 2 to column 22)",
  " (in 'rw', line 9, column 2 to column 21)",
  " (in 'rw', line 10, column 2 to column 27)",
  " (in 'rw', line 13, column 2 to column 15)",
  " (in 'rw', line 14, column 2 to column 12)",
  " (in 'rw', line 39, column 2 to column 22)",
  " (in 'rw', line 15, column 2 to column 14)",
  " (in 'rw', line 16, column 2 to column 11)",
  " (in 'rw', line 18, column 4 to column 17)",
  " (in 'rw', line 17, column 21 to line 19, column 3)",
  " (in 'rw', line 17, column 2 to line 19, column 3)",
  " (in 'rw', line 21, column 4 to column 28)",
  " (in 'rw', line 20, column 16 to line 22, column 3)",
  " (in 'rw', line 20, column 2 to line 22, column 3)",
  " (in 'rw', line 41, column 17 to column 67)",
  " (in 'rw', line 41, column 2 to column 67)",
  " (in 'rw', line 26, column 4 to column 26)",
  " (in 'rw', line 25, column 16 to line 27, column 3)",
  " (in 'rw', line 25, column 2 to line 27, column 3)",
  " (in 'rw', line 31, column 4 to column 45)",
  " (in 'rw', line 30, column 9 to line 32, column 3)",
  " (in 'rw', line 29, column 4 to column 38)",
  " (in 'rw', line 28, column 16 to line 30, column 3)",
  " (in 'rw', line 28, column 2 to line 32, column 3)",
  " (in 'rw', line 33, column 2 to column 27)",
  " (in 'rw', line 35, column 4 to column 21)",
  " (in 'rw', line 34, column 19 to line 36, column 3)",
  " (in 'rw', line 34, column 2 to line 36, column 3)",
  " (in 'rw', line 2, column 2 to column 17)",
  " (in 'rw', line 3, column 9 to column 10)",
  " (in 'rw', line 3, column 2 to column 12)",
  " (in 'rw', line 4, column 2 to column 25)",
  " (in 'rw', line 5, column 2 to column 22)",
  " (in 'rw', line 9, column 10 to column 19)",
  " (in 'rw', line 10, column 19 to column 25)",
  " (in 'rw', line 13, column 12 to column 13)",
  " (in 'rw', line 39, column 9 to column 12)"};
#include <stan_meta_header.hpp>
class model_rw final : public model_base_crtp<model_rw> {
private:
  int N;
  std::vector<double> y;
  int est_drift;
  int est_nu;
  int log_lik_1dim__;
public:
  ~model_rw() {}
  model_rw(stan::io::var_context& context__, unsigned int random_seed__ = 0,
           std::ostream* pstream__ = nullptr) : model_base_crtp(0) {
    int current_statement__ = 0;
    using local_scalar_t__ = double;
    boost::ecuyer1988 base_rng__ =
      stan::services::util::create_rng(random_seed__, 0);
    // suppress unused var warning
    (void) base_rng__;
    static constexpr const char* function__ = "model_rw_namespace::model_rw";
    // suppress unused var warning
    (void) function__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      current_statement__ = 29;
      context__.validate_dims("data initialization", "N", "int",
        std::vector<size_t>{});
      N = std::numeric_limits<int>::min();
      current_statement__ = 29;
      N = context__.vals_i("N")[(1 - 1)];
      current_statement__ = 29;
      stan::math::check_greater_or_equal(function__, "N", N, 0);
      current_statement__ = 30;
      stan::math::validate_non_negative_index("y", "N", N);
      current_statement__ = 31;
      context__.validate_dims("data initialization", "y", "double",
        std::vector<size_t>{static_cast<size_t>(N)});
      y = std::vector<double>(N, std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 31;
      y = context__.vals_r("y");
      current_statement__ = 32;
      context__.validate_dims("data initialization", "est_drift", "int",
        std::vector<size_t>{});
      est_drift = std::numeric_limits<int>::min();
      current_statement__ = 32;
      est_drift = context__.vals_i("est_drift")[(1 - 1)];
      current_statement__ = 32;
      stan::math::check_greater_or_equal(function__, "est_drift", est_drift,
        0);
      current_statement__ = 33;
      context__.validate_dims("data initialization", "est_nu", "int",
        std::vector<size_t>{});
      est_nu = std::numeric_limits<int>::min();
      current_statement__ = 33;
      est_nu = context__.vals_i("est_nu")[(1 - 1)];
      current_statement__ = 33;
      stan::math::check_greater_or_equal(function__, "est_nu", est_nu, 0);
      current_statement__ = 34;
      stan::math::validate_non_negative_index("mu", "est_drift", est_drift);
      current_statement__ = 35;
      stan::math::validate_non_negative_index("nu", "est_nu", est_nu);
      current_statement__ = 36;
      stan::math::validate_non_negative_index("pred", "N", N);
      current_statement__ = 37;
      log_lik_1dim__ = std::numeric_limits<int>::min();
      current_statement__ = 37;
      log_lik_1dim__ = (N - 1);
      current_statement__ = 37;
      stan::math::validate_non_negative_index("log_lik", "N - 1",
        log_lik_1dim__);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    num_params_r__ = 1 + est_drift + est_nu;
  }
  inline std::string model_name() const final {
    return "model_rw";
  }
  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 v2.32.2",
             "stancflags = --allow-undefined"};
  }
  template <bool propto__, bool jacobian__, typename VecR, typename VecI,
            stan::require_vector_like_t<VecR>* = nullptr,
            stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline stan::scalar_type_t<VecR>
  log_prob_impl(VecR& params_r__, VecI& params_i__, std::ostream*
                pstream__ = nullptr) const {
    using T__ = stan::scalar_type_t<VecR>;
    using local_scalar_t__ = T__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    static constexpr const char* function__ = "model_rw_namespace::log_prob";
    // suppress unused var warning
    (void) function__;
    try {
      local_scalar_t__ sigma = DUMMY_VAR__;
      current_statement__ = 1;
      sigma = in__.template read_constrain_lb<local_scalar_t__,
                jacobian__>(0, lp__);
      std::vector<local_scalar_t__> mu =
        std::vector<local_scalar_t__>(est_drift, DUMMY_VAR__);
      current_statement__ = 2;
      mu = in__.template read<std::vector<local_scalar_t__>>(est_drift);
      std::vector<local_scalar_t__> nu =
        std::vector<local_scalar_t__>(est_nu, DUMMY_VAR__);
      current_statement__ = 3;
      nu = in__.template read_constrain_lb<std::vector<local_scalar_t__>,
             jacobian__>(2, lp__, est_nu);
      std::vector<local_scalar_t__> pred =
        std::vector<local_scalar_t__>(N, DUMMY_VAR__);
      local_scalar_t__ temp = DUMMY_VAR__;
      current_statement__ = 7;
      stan::model::assign(pred, 0, "assigning variable pred",
        stan::model::index_uni(1));
      current_statement__ = 8;
      temp = 0;
      current_statement__ = 11;
      if (stan::math::logical_eq(est_drift, 1)) {
        current_statement__ = 9;
        temp = stan::model::rvalue(mu, "mu", stan::model::index_uni(1));
      }
      current_statement__ = 14;
      for (int i = 2; i <= N; ++i) {
        current_statement__ = 12;
        stan::model::assign(pred,
          (stan::model::rvalue(y, "y", stan::model::index_uni((i - 1))) +
          temp), "assigning variable pred", stan::model::index_uni(i));
      }
      {
        current_statement__ = 19;
        if (stan::math::logical_eq(est_nu, 1)) {
          current_statement__ = 17;
          lp_accum__.add(stan::math::student_t_lpdf<propto__>(nu, 3, 2, 2));
        }
        current_statement__ = 24;
        if (stan::math::logical_eq(est_nu, 0)) {
          current_statement__ = 22;
          lp_accum__.add(stan::math::normal_lpdf<propto__>(
                           stan::model::rvalue(y, "y",
                             stan::model::index_min_max(2, N)),
                           stan::model::rvalue(pred, "pred",
                             stan::model::index_min_max(2, N)), sigma));
        } else {
          current_statement__ = 20;
          lp_accum__.add(stan::math::student_t_lpdf<propto__>(
                           stan::model::rvalue(y, "y",
                             stan::model::index_min_max(2, N)), nu,
                           stan::model::rvalue(pred, "pred",
                             stan::model::index_min_max(2, N)), sigma));
        }
        current_statement__ = 25;
        lp_accum__.add(stan::math::student_t_lpdf<propto__>(sigma, 3, 0, 2));
        current_statement__ = 28;
        if (stan::math::logical_eq(est_drift, 1)) {
          current_statement__ = 26;
          lp_accum__.add(stan::math::normal_lpdf<propto__>(mu, 0, 1));
        }
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
  }
  template <typename RNG, typename VecR, typename VecI, typename VecVar,
            stan::require_vector_like_vt<std::is_floating_point,
            VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral,
            VecI>* = nullptr, stan::require_vector_vt<std::is_floating_point,
            VecVar>* = nullptr>
  inline void
  write_array_impl(RNG& base_rng__, VecR& params_r__, VecI& params_i__,
                   VecVar& vars__, const bool
                   emit_transformed_parameters__ = true, const bool
                   emit_generated_quantities__ = true, std::ostream*
                   pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    static constexpr bool propto__ = true;
    // suppress unused var warning
    (void) propto__;
    double lp__ = 0.0;
    // suppress unused var warning
    (void) lp__;
    int current_statement__ = 0;
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    constexpr bool jacobian__ = false;
    static constexpr const char* function__ =
      "model_rw_namespace::write_array";
    // suppress unused var warning
    (void) function__;
    try {
      double sigma = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 1;
      sigma = in__.template read_constrain_lb<local_scalar_t__,
                jacobian__>(0, lp__);
      std::vector<double> mu =
        std::vector<double>(est_drift,
          std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 2;
      mu = in__.template read<std::vector<local_scalar_t__>>(est_drift);
      std::vector<double> nu =
        std::vector<double>(est_nu, std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 3;
      nu = in__.template read_constrain_lb<std::vector<local_scalar_t__>,
             jacobian__>(2, lp__, est_nu);
      std::vector<double> pred =
        std::vector<double>(N, std::numeric_limits<double>::quiet_NaN());
      double temp = std::numeric_limits<double>::quiet_NaN();
      out__.write(sigma);
      out__.write(mu);
      out__.write(nu);
      if (stan::math::logical_negation(
            (stan::math::primitive_value(emit_transformed_parameters__) ||
            stan::math::primitive_value(emit_generated_quantities__)))) {
        return ;
      }
      current_statement__ = 7;
      stan::model::assign(pred, 0, "assigning variable pred",
        stan::model::index_uni(1));
      current_statement__ = 8;
      temp = 0;
      current_statement__ = 11;
      if (stan::math::logical_eq(est_drift, 1)) {
        current_statement__ = 9;
        temp = stan::model::rvalue(mu, "mu", stan::model::index_uni(1));
      }
      current_statement__ = 14;
      for (int i = 2; i <= N; ++i) {
        current_statement__ = 12;
        stan::model::assign(pred,
          (stan::model::rvalue(y, "y", stan::model::index_uni((i - 1))) +
          temp), "assigning variable pred", stan::model::index_uni(i));
      }
      if (emit_transformed_parameters__) {
        out__.write(pred);
        out__.write(temp);
      }
      if (stan::math::logical_negation(emit_generated_quantities__)) {
        return ;
      }
      Eigen::Matrix<double,-1,1> log_lik =
        Eigen::Matrix<double,-1,1>::Constant(log_lik_1dim__,
          std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 16;
      for (int n = 2; n <= N; ++n) {
        current_statement__ = 15;
        stan::model::assign(log_lik,
          stan::math::normal_lpdf<false>(
            stan::model::rvalue(y, "y", stan::model::index_uni(n)),
            stan::model::rvalue(pred, "pred", stan::model::index_uni(n)),
            sigma), "assigning variable log_lik",
          stan::model::index_uni((n - 1)));
      }
      out__.write(log_lik);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  template <typename VecVar, typename VecI,
            stan::require_vector_t<VecVar>* = nullptr,
            stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline void
  unconstrain_array_impl(const VecVar& params_r__, const VecI& params_i__,
                         VecVar& vars__, std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      local_scalar_t__ sigma = DUMMY_VAR__;
      current_statement__ = 1;
      sigma = in__.read<local_scalar_t__>();
      out__.write_free_lb(0, sigma);
      std::vector<local_scalar_t__> mu =
        std::vector<local_scalar_t__>(est_drift, DUMMY_VAR__);
      current_statement__ = 2;
      stan::model::assign(mu,
        in__.read<std::vector<local_scalar_t__>>(est_drift),
        "assigning variable mu");
      out__.write(mu);
      std::vector<local_scalar_t__> nu =
        std::vector<local_scalar_t__>(est_nu, DUMMY_VAR__);
      current_statement__ = 3;
      stan::model::assign(nu,
        in__.read<std::vector<local_scalar_t__>>(est_nu),
        "assigning variable nu");
      out__.write_free_lb(2, nu);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  template <typename VecVar, stan::require_vector_t<VecVar>* = nullptr>
  inline void
  transform_inits_impl(const stan::io::var_context& context__, VecVar&
                       vars__, std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::serializer<local_scalar_t__> out__(vars__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      current_statement__ = 1;
      context__.validate_dims("parameter initialization", "sigma", "double",
        std::vector<size_t>{});
      current_statement__ = 2;
      context__.validate_dims("parameter initialization", "mu", "double",
        std::vector<size_t>{static_cast<size_t>(est_drift)});
      current_statement__ = 3;
      context__.validate_dims("parameter initialization", "nu", "double",
        std::vector<size_t>{static_cast<size_t>(est_nu)});
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      local_scalar_t__ sigma = DUMMY_VAR__;
      current_statement__ = 1;
      sigma = context__.vals_r("sigma")[(1 - 1)];
      out__.write_free_lb(0, sigma);
      std::vector<local_scalar_t__> mu =
        std::vector<local_scalar_t__>(est_drift, DUMMY_VAR__);
      current_statement__ = 2;
      mu = context__.vals_r("mu");
      out__.write(mu);
      std::vector<local_scalar_t__> nu =
        std::vector<local_scalar_t__>(est_nu, DUMMY_VAR__);
      current_statement__ = 3;
      nu = context__.vals_r("nu");
      out__.write_free_lb(2, nu);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  inline void
  get_param_names(std::vector<std::string>& names__, const bool
                  emit_transformed_parameters__ = true, const bool
                  emit_generated_quantities__ = true) const {
    names__ = std::vector<std::string>{"sigma", "mu", "nu"};
    if (emit_transformed_parameters__) {
      std::vector<std::string> temp{"pred", "temp"};
      names__.reserve(names__.size() + temp.size());
      names__.insert(names__.end(), temp.begin(), temp.end());
    }
    if (emit_generated_quantities__) {
      std::vector<std::string> temp{"log_lik"};
      names__.reserve(names__.size() + temp.size());
      names__.insert(names__.end(), temp.begin(), temp.end());
    }
  }
  inline void
  get_dims(std::vector<std::vector<size_t>>& dimss__, const bool
           emit_transformed_parameters__ = true, const bool
           emit_generated_quantities__ = true) const {
    dimss__ = std::vector<std::vector<size_t>>{std::vector<size_t>{},
                std::vector<size_t>{static_cast<size_t>(est_drift)},
                std::vector<size_t>{static_cast<size_t>(est_nu)}};
    if (emit_transformed_parameters__) {
      std::vector<std::vector<size_t>>
        temp{std::vector<size_t>{static_cast<size_t>(N)},
             std::vector<size_t>{}};
      dimss__.reserve(dimss__.size() + temp.size());
      dimss__.insert(dimss__.end(), temp.begin(), temp.end());
    }
    if (emit_generated_quantities__) {
      std::vector<std::vector<size_t>>
        temp{std::vector<size_t>{static_cast<size_t>(log_lik_1dim__)}};
      dimss__.reserve(dimss__.size() + temp.size());
      dimss__.insert(dimss__.end(), temp.begin(), temp.end());
    }
  }
  inline void
  constrained_param_names(std::vector<std::string>& param_names__, bool
                          emit_transformed_parameters__ = true, bool
                          emit_generated_quantities__ = true) const final {
    param_names__.emplace_back(std::string() + "sigma");
    for (int sym1__ = 1; sym1__ <= est_drift; ++sym1__) {
      param_names__.emplace_back(std::string() + "mu" + '.' +
        std::to_string(sym1__));
    }
    for (int sym1__ = 1; sym1__ <= est_nu; ++sym1__) {
      param_names__.emplace_back(std::string() + "nu" + '.' +
        std::to_string(sym1__));
    }
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
        param_names__.emplace_back(std::string() + "pred" + '.' +
          std::to_string(sym1__));
      }
      param_names__.emplace_back(std::string() + "temp");
    }
    if (emit_generated_quantities__) {
      for (int sym1__ = 1; sym1__ <= log_lik_1dim__; ++sym1__) {
        param_names__.emplace_back(std::string() + "log_lik" + '.' +
          std::to_string(sym1__));
      }
    }
  }
  inline void
  unconstrained_param_names(std::vector<std::string>& param_names__, bool
                            emit_transformed_parameters__ = true, bool
                            emit_generated_quantities__ = true) const final {
    param_names__.emplace_back(std::string() + "sigma");
    for (int sym1__ = 1; sym1__ <= est_drift; ++sym1__) {
      param_names__.emplace_back(std::string() + "mu" + '.' +
        std::to_string(sym1__));
    }
    for (int sym1__ = 1; sym1__ <= est_nu; ++sym1__) {
      param_names__.emplace_back(std::string() + "nu" + '.' +
        std::to_string(sym1__));
    }
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
        param_names__.emplace_back(std::string() + "pred" + '.' +
          std::to_string(sym1__));
      }
      param_names__.emplace_back(std::string() + "temp");
    }
    if (emit_generated_quantities__) {
      for (int sym1__ = 1; sym1__ <= log_lik_1dim__; ++sym1__) {
        param_names__.emplace_back(std::string() + "log_lik" + '.' +
          std::to_string(sym1__));
      }
    }
  }
  inline std::string get_constrained_sizedtypes() const {
    return std::string("[{\"name\":\"sigma\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"mu\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(est_drift) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"parameters\"},{\"name\":\"nu\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(est_nu) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"parameters\"},{\"name\":\"pred\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(N) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},{\"name\":\"temp\",\"type\":{\"name\":\"real\"},\"block\":\"transformed_parameters\"},{\"name\":\"log_lik\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(log_lik_1dim__) + "},\"block\":\"generated_quantities\"}]");
  }
  inline std::string get_unconstrained_sizedtypes() const {
    return std::string("[{\"name\":\"sigma\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"mu\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(est_drift) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"parameters\"},{\"name\":\"nu\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(est_nu) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"parameters\"},{\"name\":\"pred\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(N) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},{\"name\":\"temp\",\"type\":{\"name\":\"real\"},\"block\":\"transformed_parameters\"},{\"name\":\"log_lik\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(log_lik_1dim__) + "},\"block\":\"generated_quantities\"}]");
  }
  // Begin method overload boilerplate
  template <typename RNG> inline void
  write_array(RNG& base_rng, Eigen::Matrix<double,-1,1>& params_r,
              Eigen::Matrix<double,-1,1>& vars, const bool
              emit_transformed_parameters = true, const bool
              emit_generated_quantities = true, std::ostream*
              pstream = nullptr) const {
    const size_t num_params__ = ((1 + est_drift) + est_nu);
    const size_t num_transformed = emit_transformed_parameters * ((N + 1));
    const size_t num_gen_quantities = emit_generated_quantities *
      (log_lik_1dim__);
    const size_t num_to_write = num_params__ + num_transformed +
      num_gen_quantities;
    std::vector<int> params_i;
    vars = Eigen::Matrix<double,-1,1>::Constant(num_to_write,
             std::numeric_limits<double>::quiet_NaN());
    write_array_impl(base_rng, params_r, params_i, vars,
      emit_transformed_parameters, emit_generated_quantities, pstream);
  }
  template <typename RNG> inline void
  write_array(RNG& base_rng, std::vector<double>& params_r, std::vector<int>&
              params_i, std::vector<double>& vars, bool
              emit_transformed_parameters = true, bool
              emit_generated_quantities = true, std::ostream*
              pstream = nullptr) const {
    const size_t num_params__ = ((1 + est_drift) + est_nu);
    const size_t num_transformed = emit_transformed_parameters * ((N + 1));
    const size_t num_gen_quantities = emit_generated_quantities *
      (log_lik_1dim__);
    const size_t num_to_write = num_params__ + num_transformed +
      num_gen_quantities;
    vars = std::vector<double>(num_to_write,
             std::numeric_limits<double>::quiet_NaN());
    write_array_impl(base_rng, params_r, params_i, vars,
      emit_transformed_parameters, emit_generated_quantities, pstream);
  }
  template <bool propto__, bool jacobian__, typename T_> inline T_
  log_prob(Eigen::Matrix<T_,-1,1>& params_r, std::ostream* pstream = nullptr) const {
    Eigen::Matrix<int,-1,1> params_i;
    return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
  }
  template <bool propto__, bool jacobian__, typename T_> inline T_
  log_prob(std::vector<T_>& params_r, std::vector<int>& params_i,
           std::ostream* pstream = nullptr) const {
    return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
  }
  inline void
  transform_inits(const stan::io::var_context& context,
                  Eigen::Matrix<double,-1,1>& params_r, std::ostream*
                  pstream = nullptr) const final {
    std::vector<double> params_r_vec(params_r.size());
    std::vector<int> params_i;
    transform_inits(context, params_i, params_r_vec, pstream);
    params_r = Eigen::Map<Eigen::Matrix<double,-1,1>>(params_r_vec.data(),
                 params_r_vec.size());
  }
  inline void
  transform_inits(const stan::io::var_context& context, std::vector<int>&
                  params_i, std::vector<double>& vars, std::ostream*
                  pstream__ = nullptr) const {
    vars.resize(num_params_r__);
    transform_inits_impl(context, vars, pstream__);
  }
  inline void
  unconstrain_array(const std::vector<double>& params_constrained,
                    std::vector<double>& params_unconstrained, std::ostream*
                    pstream = nullptr) const {
    const std::vector<int> params_i;
    params_unconstrained = std::vector<double>(num_params_r__,
                             std::numeric_limits<double>::quiet_NaN());
    unconstrain_array_impl(params_constrained, params_i,
      params_unconstrained, pstream);
  }
  inline void
  unconstrain_array(const Eigen::Matrix<double,-1,1>& params_constrained,
                    Eigen::Matrix<double,-1,1>& params_unconstrained,
                    std::ostream* pstream = nullptr) const {
    const std::vector<int> params_i;
    params_unconstrained = Eigen::Matrix<double,-1,1>::Constant(num_params_r__,
                             std::numeric_limits<double>::quiet_NaN());
    unconstrain_array_impl(params_constrained, params_i,
      params_unconstrained, pstream);
  }
};
}
using stan_model = model_rw_namespace::model_rw;
#ifndef USING_R
// Boilerplate
stan::model::model_base&
new_model(stan::io::var_context& data_context, unsigned int seed,
          std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
stan::math::profile_map& get_stan_profile_data() {
  return model_rw_namespace::profiles__;
}
#endif
#endif
