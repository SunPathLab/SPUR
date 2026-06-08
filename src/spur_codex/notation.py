from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, order=True)
class CopyNumberState:
    total_cn: int
    minor_cn: int

    def __post_init__(self) -> None:
        if self.total_cn < 0:
            raise ValueError("Total copy number must be nonnegative.")
        if self.minor_cn < 0:
            raise ValueError("Minor copy number must be nonnegative.")
        if self.minor_cn > self.total_cn:
            raise ValueError("Minor copy number cannot exceed total copy number.")
        if self.minor_cn > self.major_cn:
            raise ValueError("Minor copy number cannot exceed major copy number.")

    @property
    def major_cn(self) -> int:
        return self.total_cn - self.minor_cn

    @property
    def cn_tuple(self) -> tuple[int, int]:
        return self.total_cn, self.minor_cn

    @property
    def major_minor_tuple(self) -> tuple[int, int]:
        return self.major_cn, self.minor_cn

    @property
    def label(self) -> str:
        return f"{self.total_cn}:{self.minor_cn}"

    @classmethod
    def from_label(cls, value: str) -> "CopyNumberState":
        total_cn, minor_cn = (int(part) for part in value.split(":"))
        return cls(total_cn=total_cn, minor_cn=minor_cn)

    @classmethod
    def from_tuple(cls, cn: tuple[int, int]) -> "CopyNumberState":
        total_cn, minor_cn = cn
        return cls(total_cn=total_cn, minor_cn=minor_cn)


@dataclass(frozen=True)
class CopyNumberPair:
    sample1: CopyNumberState
    sample2: CopyNumberState

    @property
    def sample1_cn_tuple(self) -> tuple[int, int]:
        return self.sample1.cn_tuple

    @property
    def sample2_cn_tuple(self) -> tuple[int, int]:
        return self.sample2.cn_tuple

    @property
    def label(self) -> str:
        return f"{self.sample1.label}_{self.sample2.label}"

    def canonical(self) -> "CopyNumberPair":
        sample1, sample2 = sorted(
            (self.sample1, self.sample2),
            key=lambda state: (state.total_cn, state.minor_cn, state.major_cn),
            reverse=True,
        )
        return CopyNumberPair(sample1=sample1, sample2=sample2)

    @classmethod
    def from_label(cls, value: str) -> "CopyNumberPair":
        sample1_label, sample2_label = value.split("_")
        return cls(
            sample1=CopyNumberState.from_label(sample1_label),
            sample2=CopyNumberState.from_label(sample2_label),
        )


def total_minor_label(cn: tuple[int, int]) -> str:
    return CopyNumberState.from_tuple(cn).label


def canonical_cn_pair(sample1_cn: tuple[int, int], sample2_cn: tuple[int, int]) -> tuple[tuple[int, int], tuple[int, int]]:
    pair = CopyNumberPair(
        sample1=CopyNumberState.from_tuple(sample1_cn),
        sample2=CopyNumberState.from_tuple(sample2_cn),
    ).canonical()
    return pair.sample1_cn_tuple, pair.sample2_cn_tuple


def parse_between_sample_cn(value: str) -> tuple[tuple[int, int], tuple[int, int]]:
    pair = CopyNumberPair.from_label(value).canonical()
    return pair.sample1_cn_tuple, pair.sample2_cn_tuple


def truncal_key(total_cn: int, minor_cn: int, tree_shape: str) -> str:
    return f"truncal_{total_cn}:{minor_cn}_{tree_shape}"


def to_major_minor_tuple(cn: tuple[int, int]) -> tuple[int, int]:
    return CopyNumberState.from_tuple(cn).major_minor_tuple
