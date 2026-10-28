#include "var_op_function.h"
#include <cmath>

// external initialization required
std::unique_ptr<fft_op> fft_obj;

/**
 * @brief Computes L2 error metrics between two real-valued vectors.
 *
 * Calculates:
 * - RMSE (root mean square error)
 * - Relative L2 error (normalized by signal energy)
 *
 * @param orig Original signal vector.
 * @param rec  Reconstructed signal vector.
 * @return Pair {RMSE, relative L2 error}.
 *
 * @throws std::runtime_error if vector sizes differ.
 */
std::pair<double, double> compute_l2_error(const std::vector<double>& orig, const std::vector<double>& rec) {
    if (orig.size() != rec.size()) {
        throw std::runtime_error("Two vectors must have same size in compute_l2_error!");
    }
    double sum_sq_diff = 0.0;
    double sum_sq_orig = 0.0;
    size_t n = orig.size();
    if (n == 0) return {0.0, 0.0};
    for (size_t i = 0; i < n; ++i) {
        double diff = orig[i] - rec[i];
        sum_sq_diff += diff * diff;
        sum_sq_orig += orig[i] * orig[i];
    }
    double rmse = std::sqrt(sum_sq_diff / n);
    double rel_err = (sum_sq_orig > 0.0) ? std::sqrt(sum_sq_diff / sum_sq_orig) : 0.0;
    return std::make_pair(rmse, rel_err);
}

/**
 * @brief Computes L2 error metrics between two complex-valued vectors.
 *
 * Calculates:
 * - RMSE based on complex magnitude difference
 * - Relative L2 error normalized by original signal energy
 *
 * @param orig Original complex signal vector.
 * @param rec  Reconstructed complex signal vector.
 * @return Pair {RMSE, relative L2 error}.
 *
 * @throws std::runtime_error if vector sizes differ.
 */
std::pair<double, double> compute_l2_error(const std::vector<std::complex<double>>& orig, const std::vector<std::complex<double>>& rec) {
    if (orig.size() != rec.size()) {
        throw std::runtime_error("Two vectors must have same size in compute_l2_error!");
    }
    
    double sum_sq_diff = 0.0;
    double sum_sq_orig = 0.0;
    size_t n = orig.size();
    if (n == 0) return {0.0, 0.0};
    for (size_t i = 0; i < n; ++i) {
        double diff = std::abs(orig[i] - rec[i]);
        sum_sq_diff += diff * diff;
        sum_sq_orig += std::norm(orig[i]);
    }
    double rmse = std::sqrt(sum_sq_diff / n);
    double rel_err = (sum_sq_orig > 0.0) ? std::sqrt(sum_sq_diff / sum_sq_orig) : 0.0;
    return std::make_pair(rmse, rel_err);
}

// ---------------------------------------------------------------------
// Implementation 1D fft tools.
// ---------------------------------------------------------------------

//class fft1d::constructor
fft1d::fft1d(size_t N)
{
    if (N == 0)
    {
        throw std::invalid_argument(
            "fft1d: N must be positive."
        );
    }
    N_ = N;
    in_buffer_.resize(N_);
    out_buffer_.resize(N_);

    forward_plan_ =
        fftw_plan_dft_1d(
            N_,
            reinterpret_cast<fftw_complex*>(
                in_buffer_.data()),
            reinterpret_cast<fftw_complex*>(
                out_buffer_.data()),
            FFTW_FORWARD,
            FFTW_MEASURE
        );

    if (!forward_plan_)
    {
        throw std::runtime_error(
            "fft1d: failed to create forward plan."
        );
    }

    backward_plan_ =
        fftw_plan_dft_1d(
            N_,
            reinterpret_cast<fftw_complex*>(
                in_buffer_.data()),
            reinterpret_cast<fftw_complex*>(
                out_buffer_.data()),
            FFTW_BACKWARD,
            FFTW_MEASURE
        );

    if (!backward_plan_)
    {
        fftw_destroy_plan(
            forward_plan_);

        throw std::runtime_error(
            "fft1d: failed to create inverse plan."
        );
    }
}
//class fft1d::destructor
fft1d::~fft1d()
{
    fftw_destroy_plan(
        forward_plan_);

    fftw_destroy_plan(
        backward_plan_);
}
// forward transform
const std::vector<std::complex<double>>&
fft1d::fft(
    const std::vector<std::complex<double>>& input)
{
    if (input.size() != N_)
    {
        throw std::invalid_argument(
            "fft1d::fft: wrong input size."
        );
    }

    std::copy(
        input.begin(),
        input.end(),
        in_buffer_.begin()
    );

    fftw_execute(
        forward_plan_);

    return out_buffer_;
}
// backward transform
const std::vector<std::complex<double>>&
fft1d::ifft(
    const std::vector<std::complex<double>>& input)
{
    if (input.size() != N_)
    {
        throw std::invalid_argument(
            "fft1d::ifft: wrong input size."
        );
    }

    std::copy(
        input.begin(),
        input.end(),
        in_buffer_.begin()
    );

    fftw_execute(
        backward_plan_);

    const double scale =
        1.0 / static_cast<double>(N_);

    for (auto& z : out_buffer_)
    {
        z *= scale;
    }

    return out_buffer_;
}

// ---------------------------------------------------------------------
// Implementation 2D fft tools.
// ---------------------------------------------------------------------

fft2d::fft2d(size_t Nx, size_t Ny)
    : Nx_(Nx), Ny_(Ny), total_(Nx * Ny)
{
    if (Nx == 0 || Ny == 0) {
        throw std::invalid_argument("fft2d: Nx and Ny must be positive.");
    }

    in_buffer_.resize(total_);
    out_buffer_.resize(total_);

    // Create forward plan (FFTW_FORWARD)
    int n[2] = { static_cast<int>(Nx_), static_cast<int>(Ny_) };
    forward_plan_ = fftw_plan_dft_2d(
        n[0], n[1],
        reinterpret_cast<fftw_complex*>(in_buffer_.data()),
        reinterpret_cast<fftw_complex*>(out_buffer_.data()),
        FFTW_FORWARD,
        FFTW_MEASURE
    );
    if (!forward_plan_) {
        throw std::runtime_error("fft2d: failed to create forward plan.");
    }

    // Create backward plan (FFTW_BACKWARD)
    backward_plan_ = fftw_plan_dft_2d(
        n[0], n[1],
        reinterpret_cast<fftw_complex*>(in_buffer_.data()),
        reinterpret_cast<fftw_complex*>(out_buffer_.data()),
        FFTW_BACKWARD,
        FFTW_MEASURE
    );
    if (!backward_plan_) {
        fftw_destroy_plan(forward_plan_);
        throw std::runtime_error("fft2d: failed to create inverse plan.");
    }
}

fft2d::~fft2d()
{
    fftw_destroy_plan(forward_plan_);
    fftw_destroy_plan(backward_plan_);
}

const std::vector<std::complex<double>>&
fft2d::fft2(const std::vector<std::complex<double>>& input)
{
    if (input.size() != total_) {
        throw std::invalid_argument("fft2d::fft2: wrong input size.");
    }

    std::copy(input.begin(), input.end(), in_buffer_.begin());
    fftw_execute(forward_plan_);
    return out_buffer_;
}

const std::vector<std::complex<double>>&
fft2d::ifft2(const std::vector<std::complex<double>>& input)
{
    if (input.size() != total_) {
        throw std::invalid_argument("fft2d::ifft2: wrong input size.");
    }

    std::copy(input.begin(), input.end(), in_buffer_.begin());
    fftw_execute(backward_plan_);

    const double scale = 1.0 / static_cast<double>(total_);
    for (auto& z : out_buffer_) {
        z *= scale;
    }
    return out_buffer_;
}

