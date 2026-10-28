"""Run CoRAL's current CLI with discordant-arc strokes at quarter width."""

from coral import plot_amplicons
from coral.cli import coral_app


ARC_LINEWIDTH_SCALE = 0.25
_original_linewidth = plot_amplicons.get_discordant_edge_linewidth


def _thin_discordant_edge_linewidth(
    edge_read_count: float,
    max_read_count: float,
) -> float:
    return ARC_LINEWIDTH_SCALE * _original_linewidth(
        edge_read_count,
        max_read_count,
    )


plot_amplicons.get_discordant_edge_linewidth = (
    _thin_discordant_edge_linewidth
)


if __name__ == "__main__":
    coral_app()
