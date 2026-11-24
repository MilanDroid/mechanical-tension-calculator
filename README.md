# Mechanical tension calculator

*Cálculo mecánico de conductores*

A browser tool that calculates the tension, safety factor and sag of an overhead power line conductor under each standard load case. I built it as a university project at UNAH (Universidad Nacional Autónoma de Honduras).

**Try it:** [milandroid.dev/mechanical-tension-calculator](https://milandroid.dev/mechanical-tension-calculator/). The interface is in Spanish.

## What you enter

- **Conductor:** diameter (mm), cross-section (mm²), weight (kg/m), breaking load (kg), modulus of elasticity (kg/mm²) and thermal expansion coefficient.
- **Line:** altitude (m above sea level) and average span.
- **Cable type:** phase conductor or earth wire.

## What it calculates

**1. The starting condition, chosen by altitude.**

| Zone | Altitude | Condition | Extra load per metre |
|---|---|---|---|
| A | below 500 m | −5 °C with wind | 60·d/1000 kg for d ≤ 16 mm, 50·d/1000 kg above that |
| B | 500 m to under 1000 m | −15 °C with ice | 0.18·√d kg |
| C | 1000 m and above | −20 °C with ice | 0.36·√d kg |

*d* is the conductor diameter in mm. In the starting condition, the cable carries its maximum allowed tension: the breaking load divided by 3 for a phase conductor, or by 5 for an earth wire.

**2. Every other load case, from the change-of-state equation.** For each case, the calculator solves the equation (a cubic) for the new tension:

| Case | Condition |
|---|---|
| D | 15 °C with wind |
| E | 50 °C, no extra load |
| F | 0 °C, with ice in zones B and C (no extra load in zone A) |
| G | 15 °C, no extra load |
| H | −5 °C, no extra load |
| I | The starting temperature (−5, −15 or −20 °C), no extra load |
| J | −5 °C with wind |

**3. The results.** For every case, it shows the tension (kg), the safety factor, and the vertical and inclined sag (m). It also reports two checks as a percentage of the breaking load:

- **TCD** (*tensión de cada día*, everyday tension): the tension in case G.
- **THF** (*tensión de horas frías*, cold-hours tension): the tension in case H.

## Run it locally

There's no build step. Clone the repository and open `index.html` in a browser.

| File | Contents |
|---|---|
| `index.html` | The form and the results table |
| `resources/index.js` | The calculations, including the cubic solver |
| `resources/styles.css` | Styles |
