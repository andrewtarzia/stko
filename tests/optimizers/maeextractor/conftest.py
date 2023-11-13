import pytest

from .case_data import CaseData

_energies = [
    (1649.50158691406, 1),
    (1722.63293457031, 2),
    (1607.67626953125, 3),
    (1604.86291503906, 4),
    (1588.51037597656, 5),
    (1577.38671875, 6),
    (1628.57958984375, 7),
    (1579.61926269531, 8),
    (1621.49975585938, 9),
    (1627.87585449219, 10),
    (1615.9462890625, 11),
    (1600.01208496094, 12),
    (1601.82263183594, 13),
    (1669.08312988281, 14),
    (1677.03796386719, 15),
    (1598.83813476562, 16),
    (1623.77331542969, 17),
    (1601.51000976562, 18),
    (1595.39501953125, 19),
    (1601.62170410156, 20),
    (1585.79992675781, 21),
    (1692.82690429688, 22),
    (1641.13891601562, 23),
    (1648.33898925781, 24),
    (1701.05114746094, 25),
    (1609.71362304688, 26),
    (1668.08581542969, 27),
    (1642.79821777344, 28),
    (1626.81860351562, 29),
    (1598.40185546875, 30),
    (1596.69543457031, 31),
    (1670.2080078125, 32),
    (1608.22106933594, 33),
    (1657.68298339844, 34),
    (1616.05322265625, 35),
    (1574.47778320312, 36),
    (1578.36291503906, 37),
    (1610.46337890625, 38),
    (1648.68884277344, 39),
    (1607.71716308594, 40),
    (1649.50402832031, 41),
    (1672.12194824219, 42),
    (1610.67810058594, 43),
    (1630.57885742188, 44),
    (1602.98059082031, 45),
    (1638.20202636719, 46),
    (1667.32836914062, 47),
    (1581.57482910156, 48),
    (1606.1328125, 49),
    (1649.66223144531, 50),
    (1610.92443847656, 51),
    (1609.87719726562, 52),
    (1650.9794921875, 53),
    (1668.99499511719, 54),
    (1700.49877929688, 55),
    (1674.94970703125, 56),
    (1642.20751953125, 57),
    (1702.6923828125, 58),
    (1672.70031738281, 59),
    (1672.16320800781, 60),
    (1670.28076171875, 61),
    (1686.03771972656, 62),
    (1664.91577148438, 63),
    (1647.33508300781, 64),
    (1714.8671875, 65),
    (1696.59802246094, 66),
    (1645.01538085938, 67),
    (1597.09240722656, 68),
    (1613.04235839844, 69),
    (1622.42114257812, 70),
    (1605.45751953125, 71),
    (1602.25036621094, 72),
    (1594.09802246094, 73),
    (1654.03540039062, 74),
    (1609.90295410156, 75),
    (1563.08666992188, 76),
    (1652.13793945312, 77),
    (1565.83251953125, 78),
    (1640.81762695312, 79),
    (1607.97631835938, 80),
    (1649.76391601562, 81),
    (1665.44958496094, 82),
    (1651.71862792969, 83),
    (1687.31689453125, 84),
    (1679.642578125, 85),
    (1676.05688476562, 86),
    (1660.06896972656, 87),
    (1610.8662109375, 88),
    (1718.63208007812, 89),
    (1694.84851074219, 90),
    (1719.84423828125, 91),
    (1698.97583007812, 92),
    (1683.58813476562, 93),
    (1661.88049316406, 94),
    (1700.79724121094, 95),
    (1682.1865234375, 96),
    (1639.63977050781, 97),
    (1629.56079101562, 98),
    (1623.24731445312, 99),
    (1711.85949707031, 100),
    (1711.732421875, 101),
]


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            run_name="test",
            n=1,
            mae_path="test-out.mae",
            energies=_energies,
            min_energy=1563.08666992188,
            path="test-outEXTRACTED_76.mae",
            name=name,
        ),
        lambda name: CaseData(
            run_name="test",
            n=4,
            mae_path="test-out.mae",
            energies=_energies,
            min_energy=1563.08666992188,
            path="test-outEXTRACTED_76_conf_0.mae",
            name=name,
        ),
    ),
)
def case_data(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )