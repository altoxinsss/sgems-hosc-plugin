from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext
from math import comb

import numpy as np


LOCATION = 0.25
SCALE = 0.05
HIGH_THRESHOLD = 0.25
ANALYTIC_DECIMAL_PRECISION = 80
ANALYTIC_DECIMAL_TOLERANCE = Decimal("1e-45")


def _as_decimal(value: object) -> Decimal:
    """Convert a frozen analytic value without going through binary formatting."""
    if isinstance(value, Decimal):
        return value
    return Decimal(str(value))


def _decimal_ratio(numerator: int, denominator: int) -> Decimal:
    with localcontext() as context:
        context.prec = ANALYTIC_DECIMAL_PRECISION
        return +(Decimal(numerator) / Decimal(denominator))


def _decimal_cumulants(
    values: tuple[Decimal, ...], probabilities: tuple[Decimal, ...], max_order: int
) -> dict[int, Decimal]:
    """Evaluate population cumulants with a high-precision independent recurrence."""
    with localcontext() as context:
        context.prec = ANALYTIC_DECIMAL_PRECISION
        moments = {
            order: sum(
                probability * (Decimal(1) if order == 0 else value**order)
                for value, probability in zip(values, probabilities, strict=True)
            )
            for order in range(max_order + 1)
        }
        cumulants: dict[int, Decimal] = {1: moments[1]}
        for order in range(2, max_order + 1):
            cumulant = moments[order]
            for lower_order in range(1, order):
                cumulant -= Decimal(comb(order - 1, lower_order - 1)) * cumulants[lower_order] * moments[
                    order - lower_order
                ]
            cumulants[order] = +cumulant
        return cumulants


@dataclass(frozen=True)
class FiniteDiscreteDistribution:
    """An arbitrary finite distribution used by a controlled moment experiment."""

    name: str
    values: tuple[float, ...]
    probabilities: tuple[float, ...]
    # The production sampling path intentionally uses float arrays.  Frozen
    # positive controls additionally carry exact/high-precision support values
    # so their population-moment gate is independent of binary floating point.
    exact_values: tuple[Decimal, ...] | None = None
    exact_probabilities: tuple[Decimal, ...] | None = None

    def __post_init__(self) -> None:
        if not self.values or len(self.values) != len(self.probabilities):
            raise ValueError("Values and probabilities must be non-empty and aligned")
        if not np.isfinite(np.asarray(self.values, dtype=float)).all():
            raise ValueError("Distribution values must be finite")
        if any(value < 0.0 for value in self.probabilities):
            raise ValueError("Probabilities must be non-negative")
        if not np.isclose(sum(self.probabilities), 1.0, rtol=0.0, atol=1.0e-14):
            raise ValueError("Probabilities must sum to one")
        if (self.exact_values is None) != (self.exact_probabilities is None):
            raise ValueError("Exact values and probabilities must be provided together")
        if self.exact_values is not None and self.exact_probabilities is not None:
            if len(self.exact_values) != len(self.values) or len(self.exact_probabilities) != len(self.probabilities):
                raise ValueError("Exact support must align with the floating-point support")
            with localcontext() as context:
                context.prec = ANALYTIC_DECIMAL_PRECISION
                if abs(sum(self.exact_probabilities) - Decimal(1)) > ANALYTIC_DECIMAL_TOLERANCE:
                    raise ValueError("Exact probabilities must sum to one")

    def raw_moment(self, order: int) -> float:
        """Return an exact population raw moment of X before affine scaling."""
        if order < 0:
            raise ValueError("Moment order must be non-negative")
        values = np.asarray(self.values, dtype=float)
        probabilities = np.asarray(self.probabilities, dtype=float)
        return float(np.sum(probabilities * values**order))

    def cumulants(self, max_order: int = 6) -> dict[int, float]:
        """Evaluate exact population cumulants by the general moment recursion."""
        if max_order < 1:
            raise ValueError("Maximum order must be positive")
        moments = {order: self.raw_moment(order) for order in range(max_order + 1)}
        cumulants: dict[int, float] = {1: moments[1]}
        for order in range(2, max_order + 1):
            value = moments[order]
            for lower_order in range(1, order):
                value -= comb(order - 1, lower_order - 1) * cumulants[lower_order] * moments[order - lower_order]
            cumulants[order] = float(value)
        return cumulants

    def affine_cumulants(
        self,
        *,
        location: float = LOCATION,
        scale: float = SCALE,
        max_order: int = 6,
    ) -> dict[int, float]:
        """Return cumulants of Z = location + scale X through ``max_order``."""
        source = self.cumulants(max_order=max_order)
        transformed = {1: float(location + scale * source[1])}
        transformed.update({order: float(scale**order * source[order]) for order in range(2, max_order + 1)})
        return transformed

    def _analytic_support(self) -> tuple[tuple[Decimal, ...], tuple[Decimal, ...]]:
        if self.exact_values is not None and self.exact_probabilities is not None:
            return self.exact_values, self.exact_probabilities
        return (
            tuple(_as_decimal(value) for value in self.values),
            tuple(_as_decimal(probability) for probability in self.probabilities),
        )

    def high_precision_cumulants(self, max_order: int = 6) -> dict[int, Decimal]:
        """Return an independent high-precision population cumulant check."""
        if max_order < 1:
            raise ValueError("Maximum order must be positive")
        values, probabilities = self._analytic_support()
        return _decimal_cumulants(values, probabilities, max_order)

    def high_precision_affine_cumulants(
        self,
        *,
        location: float = LOCATION,
        scale: float = SCALE,
        max_order: int = 6,
    ) -> dict[int, Decimal]:
        """Scale population cumulants analytically for ``Z = location + scale X``."""
        with localcontext() as context:
            context.prec = ANALYTIC_DECIMAL_PRECISION
            source = self.high_precision_cumulants(max_order=max_order)
            location_decimal = _as_decimal(location)
            scale_decimal = _as_decimal(scale)
            transformed = {1: +(location_decimal + scale_decimal * source[1])}
            transformed.update({order: +(scale_decimal**order * source[order]) for order in range(2, max_order + 1)})
            return transformed

    def high_precision_direct_affine_cumulants(
        self,
        *,
        location: float = LOCATION,
        scale: float = SCALE,
        max_order: int = 6,
    ) -> dict[int, Decimal]:
        """Recompute transformed cumulants directly from high-precision support."""
        with localcontext() as context:
            context.prec = ANALYTIC_DECIMAL_PRECISION
            values, probabilities = self._analytic_support()
            location_decimal = _as_decimal(location)
            scale_decimal = _as_decimal(scale)
            transformed_values = tuple(+(location_decimal + scale_decimal * value) for value in values)
            return _decimal_cumulants(transformed_values, probabilities, max_order)

    def high_precision_high_probability(
        self,
        *,
        high_threshold: float = HIGH_THRESHOLD,
        location: float = LOCATION,
        scale: float = SCALE,
    ) -> Decimal:
        """Return the affine-scale high-value probability from the analytic support."""
        with localcontext() as context:
            context.prec = ANALYTIC_DECIMAL_PRECISION
            values, probabilities = self._analytic_support()
            location_decimal = _as_decimal(location)
            scale_decimal = _as_decimal(scale)
            threshold_decimal = _as_decimal(high_threshold)
            return +sum(
                probability
                for value, probability in zip(values, probabilities, strict=True)
                if location_decimal + scale_decimal * value >= threshold_decimal
            )

    def standardized_cumulant(self, order: int) -> float:
        cumulants = self.cumulants(max_order=order)
        variance = cumulants[2]
        if variance <= 0.0:
            raise ValueError("Standardized cumulants require positive variance")
        return float(cumulants[order] / variance ** (order / 2.0))

    def high_probability(
        self,
        *,
        high_threshold: float = HIGH_THRESHOLD,
        location: float = LOCATION,
        scale: float = SCALE,
    ) -> float:
        """Return the exact P(Z >= high_threshold) after the affine transform."""
        if scale == 0.0:
            raise ValueError("Affine scale must be non-zero")
        transformed = location + scale * np.asarray(self.values, dtype=float)
        probabilities = np.asarray(self.probabilities, dtype=float)
        return float(np.sum(probabilities[transformed >= high_threshold]))

    def sample_windows(
        self,
        rng: np.random.Generator,
        window_count: int,
        samples_per_window: int,
        *,
        location: float = LOCATION,
        scale: float = SCALE,
    ) -> np.ndarray:
        """Draw independent windows after applying Z = location + scale X."""
        if window_count <= 0 or samples_per_window <= 0:
            raise ValueError("Window count and sample count must be positive")
        sampled = rng.choice(
            np.asarray(self.values, dtype=float),
            size=(window_count, samples_per_window),
            p=np.asarray(self.probabilities, dtype=float),
        )
        return location + scale * sampled


def symmetric_magnitude_distribution(
    name: str,
    magnitudes: tuple[float, ...],
    probabilities: tuple[float, ...],
    *,
    exact_magnitudes: tuple[Decimal, ...] | None = None,
    exact_probabilities: tuple[Decimal, ...] | None = None,
) -> FiniteDiscreteDistribution:
    """Build a finite distribution with independently equiprobable signs."""
    if any(value <= 0.0 for value in magnitudes):
        raise ValueError("Symmetric magnitudes must be strictly positive")
    values: list[float] = []
    signed_probabilities: list[float] = []
    for magnitude, probability in zip(magnitudes, probabilities, strict=True):
        values.extend((-float(magnitude), float(magnitude)))
        signed_probabilities.extend((float(probability) / 2.0, float(probability) / 2.0))
    exact_values: tuple[Decimal, ...] | None = None
    exact_signed_probabilities: tuple[Decimal, ...] | None = None
    if (exact_magnitudes is None) != (exact_probabilities is None):
        raise ValueError("Exact magnitudes and probabilities must be provided together")
    if exact_magnitudes is not None and exact_probabilities is not None:
        if len(exact_magnitudes) != len(magnitudes) or len(exact_probabilities) != len(probabilities):
            raise ValueError("Exact magnitude support must align with floating-point support")
        exact_values = tuple(
            value
            for magnitude in exact_magnitudes
            for value in (magnitude.copy_negate(), magnitude)
        )
        with localcontext() as context:
            context.prec = ANALYTIC_DECIMAL_PRECISION
            exact_signed_probabilities = tuple(
                value
                for probability in exact_probabilities
                for value in (+(probability / Decimal(2)), +(probability / Decimal(2)))
            )
    return FiniteDiscreteDistribution(
        name,
        tuple(values),
        tuple(signed_probabilities),
        exact_values=exact_values,
        exact_probabilities=exact_signed_probabilities,
    )


@dataclass(frozen=True)
class ControlledExperiment:
    """Two classes with matched population cumulants below one target order."""

    name: str
    class_a: FiniteDiscreteDistribution
    class_b: FiniteDiscreteDistribution
    baseline_model: str
    augmented_model: str
    differing_order: int
    matched_through_order: int

    def moment_rows(self) -> list[dict[str, object]]:
        rows: list[dict[str, object]] = []
        for label, distribution in (("A", self.class_a), ("B", self.class_b)):
            source_cumulants = distribution.cumulants(max_order=6)
            transformed_cumulants = distribution.affine_cumulants(max_order=6)
            row: dict[str, object] = {
                "experiment": self.name,
                "class": label,
                "distribution": distribution.name,
                "high_probability": distribution.high_probability(),
            }
            for order in range(1, 7):
                row[f"kappa{order}_x"] = source_cumulants[order]
                row[f"kappa{order}_z"] = transformed_cumulants[order]
            for order in range(2, 7):
                row[f"lambda{order}"] = distribution.standardized_cumulant(order)
            rows.append(row)
        return rows

    def analytic_checks(self, tolerance: float = 1.0e-12) -> dict[str, object]:
        cumulants_a = self.class_a.cumulants(max_order=6)
        cumulants_b = self.class_b.cumulants(max_order=6)
        transformed_a = self.class_a.affine_cumulants(max_order=6)
        transformed_b = self.class_b.affine_cumulants(max_order=6)
        matched = all(
            np.isclose(cumulants_a[order], cumulants_b[order], rtol=0.0, atol=tolerance)
            for order in range(1, self.matched_through_order + 1)
        )
        transformed_matched = all(
            np.isclose(transformed_a[order], transformed_b[order], rtol=0.0, atol=tolerance)
            for order in range(1, self.matched_through_order + 1)
        )
        high_probability_matches = bool(
            np.isclose(self.class_a.high_probability(), self.class_b.high_probability(), rtol=0.0, atol=tolerance)
        )
        difference = not np.isclose(
            cumulants_a[self.differing_order], cumulants_b[self.differing_order], rtol=0.0, atol=tolerance
        )
        high_precision_a = self.class_a.high_precision_cumulants(max_order=6)
        high_precision_b = self.class_b.high_precision_cumulants(max_order=6)
        high_precision_transformed_a = self.class_a.high_precision_affine_cumulants(max_order=6)
        high_precision_transformed_b = self.class_b.high_precision_affine_cumulants(max_order=6)
        direct_transformed_a = self.class_a.high_precision_direct_affine_cumulants(max_order=6)
        direct_transformed_b = self.class_b.high_precision_direct_affine_cumulants(max_order=6)

        def decimal_close(left: Decimal, right: Decimal) -> bool:
            return abs(left - right) <= ANALYTIC_DECIMAL_TOLERANCE

        high_precision_matched = all(
            decimal_close(high_precision_a[order], high_precision_b[order])
            for order in range(1, self.matched_through_order + 1)
        )
        high_precision_transformed_matched = all(
            decimal_close(high_precision_transformed_a[order], high_precision_transformed_b[order])
            for order in range(1, self.matched_through_order + 1)
        )
        high_precision_difference = not decimal_close(
            high_precision_a[self.differing_order], high_precision_b[self.differing_order]
        )
        high_precision_high_probability_match = decimal_close(
            self.class_a.high_precision_high_probability(), self.class_b.high_precision_high_probability()
        )
        affine_scaling_pass = all(
            decimal_close(high_precision_transformed_a[order], direct_transformed_a[order])
            and decimal_close(high_precision_transformed_b[order], direct_transformed_b[order])
            for order in range(1, 7)
        )
        analytic_gate_pass = bool(
            matched
            and transformed_matched
            and high_probability_matches
            and difference
            and high_precision_matched
            and high_precision_transformed_matched
            and high_precision_difference
            and high_precision_high_probability_match
            and affine_scaling_pass
        )
        return {
            "experiment": self.name,
            "matched_through_order": self.matched_through_order,
            "matching_pass": bool(matched),
            "transformed_matching_pass": bool(transformed_matched),
            "high_probability_match": high_probability_matches,
            "differing_order": self.differing_order,
            "difference_pass": bool(difference),
            "high_precision_matching_pass": bool(high_precision_matched),
            "high_precision_transformed_matching_pass": bool(high_precision_transformed_matched),
            "high_precision_high_probability_match": bool(high_precision_high_probability_match),
            "high_precision_difference_pass": bool(high_precision_difference),
            "affine_scaling_pass": bool(affine_scaling_pass),
            "reported_cumulant_orders": list(range(1, 7)),
            "higher_order_cumulants_reported": bool(self.differing_order <= 6),
            "analytic_gate_pass": analytic_gate_pass,
            "class_a_cumulant_x": cumulants_a[self.differing_order],
            "class_b_cumulant_x": cumulants_b[self.differing_order],
            "class_a_standardized": self.class_a.standardized_cumulant(self.differing_order),
            "class_b_standardized": self.class_b.standardized_cumulant(self.differing_order),
        }


def controlled_experiments() -> tuple[ControlledExperiment, ControlledExperiment, ControlledExperiment]:
    """Return frozen fourth-, fifth-, and sixth-order matching-moment pairs."""
    fourth = ControlledExperiment(
        name="fourth_order_matching_moments",
        class_a=symmetric_magnitude_distribution("A_abs_1", (1.0,), (1.0,)),
        class_b=symmetric_magnitude_distribution("B_abs_0_5_or_2", (0.5, 2.0), (0.8, 0.2)),
        baseline_model="M3",
        augmented_model="M4",
        differing_order=4,
        matched_through_order=3,
    )
    fifth = ControlledExperiment(
        name="fifth_order_matching_moments",
        class_a=FiniteDiscreteDistribution(
            "A_minus2_1_3",
            (-2.0, 1.0, 3.0),
            (7.0 / 30.0, 20.0 / 30.0, 3.0 / 30.0),
            exact_values=(Decimal("-2"), Decimal("1"), Decimal("3")),
            exact_probabilities=(_decimal_ratio(7, 30), _decimal_ratio(20, 30), _decimal_ratio(3, 30)),
        ),
        class_b=FiniteDiscreteDistribution(
            "B_minus3_minus1_0_2",
            (-3.0, -1.0, 0.0, 2.0),
            (2.0 / 30.0, 5.0 / 30.0, 10.0 / 30.0, 13.0 / 30.0),
            exact_values=(Decimal("-3"), Decimal("-1"), Decimal("0"), Decimal("2")),
            exact_probabilities=(
                _decimal_ratio(2, 30),
                _decimal_ratio(5, 30),
                _decimal_ratio(10, 30),
                _decimal_ratio(13, 30),
            ),
        ),
        baseline_model="M4",
        augmented_model="M5",
        differing_order=5,
        matched_through_order=4,
    )
    with localcontext() as context:
        context.prec = ANALYTIC_DECIMAL_PRECISION
        sqrt_half = _decimal_ratio(1, 2).sqrt()
        sqrt_three = Decimal(3).sqrt()
        sqrt_two_thirds = _decimal_ratio(2, 3).sqrt()
        sqrt_four = Decimal(4).sqrt()
    sixth = ControlledExperiment(
        name="sixth_order_matching_moments",
        class_a=symmetric_magnitude_distribution(
            "A_x2_0_5_or_3",
            (float(sqrt_half), float(sqrt_three)),
            (0.8, 0.2),
            exact_magnitudes=(sqrt_half, sqrt_three),
            exact_probabilities=(Decimal("0.8"), Decimal("0.2")),
        ),
        class_b=symmetric_magnitude_distribution(
            "B_x2_2_3_or_4",
            (float(sqrt_two_thirds), float(sqrt_four)),
            (0.9, 0.1),
            exact_magnitudes=(sqrt_two_thirds, sqrt_four),
            exact_probabilities=(Decimal("0.9"), Decimal("0.1")),
        ),
        baseline_model="M5",
        augmented_model="M6",
        differing_order=6,
        matched_through_order=5,
    )
    return fourth, fifth, sixth
