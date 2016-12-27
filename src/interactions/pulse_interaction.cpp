#include "pulse_interaction.h"

PulseInteraction::PulseInteraction(const std::shared_ptr<const DotVector> &dots,
                                   const std::unique_ptr<const Pulse> pulse)
    : InteractionTable(dots), pulse(std::move(pulse))
{
}
