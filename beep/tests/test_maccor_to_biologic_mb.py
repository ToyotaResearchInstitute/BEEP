# Copyright [2020] [Toyota Research Institute]
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Unit tests for maccor protcol files to biologic modulo bat protcol files"""

import os
import unittest
import xmltodict
from collections import OrderedDict
from monty.tempfile import ScratchDir

from beep.protocol.maccor_to_biologic_mb import MaccorToBiologicMb, convert_diagnostic_v5_multi_techniques

TEST_DIR = os.path.dirname(__file__)
TEST_FILE_DIR = os.path.join(TEST_DIR, "test_files")


class ConversionTest(unittest.TestCase):
    maxDiff = None

    def maccor_values_to_biologic_value_and_unit_test(self, func, tests):
        for value_str, expected_value_str, expected_unit in tests:
            actual_value, actual_unit = func(value_str)            
            self.assertEqual(actual_value, expected_value_str)
            self.assertEqual(actual_unit, expected_unit)

    def test_convert_volts(self):
        converter = MaccorToBiologicMb()
        tests = [
            ("0.1429", "142.900", "mV"),
            ("0.1429e3", "142.900", "V"),
            ("159.3624", "159362.400", "mV"),
            ("152.9", "152.900",  "V")
        ]
        self.maccor_values_to_biologic_value_and_unit_test(
            converter._convert_volts,
            tests,
        )
    
    def test_convert_amps(self):
        converter = MaccorToBiologicMb()
        tests = [
            ("0.1429", "142.900", "mA"),
            ("1.23", "1.230", "A"),
            ("152.9", "152.900",  "A"),
            ("1.2e-4", "120.000", "\N{Micro Sign}A")
        ]
        self.maccor_values_to_biologic_value_and_unit_test(
            converter._convert_amps,
            tests,
        )

    def test_convert_watts(self):
        converter = MaccorToBiologicMb()
        tests = [
            ("0.1429", "142.900", "mW"),
            ("1.23", "1.230", "W"),
            ("152.9", "152.900",  "W"),
            ("1.2e-5", "12.000", "\N{Micro Sign}W")
        ]
        self.maccor_values_to_biologic_value_and_unit_test(
            converter._convert_watts,
            tests,
        )

    def test_convert_ohms(self):
        converter = MaccorToBiologicMb()
        tests = [
            ("0.1429", "142.900", "mOhms"),
            ("1.459e4", "14.590", "kOhms"),
            ("152.9", "152.900",  "Ohms"),
            ("1.2e-4", "120.000", "\N{Micro Sign}Ohms")
        ]
        self.maccor_values_to_biologic_value_and_unit_test(
            converter._convert_ohms,
            tests,
        )
    
    def test_convert_time(self):
        converter = MaccorToBiologicMb()
        tests = [
            ("::.01", "10.000", "ms"),
            ("03::", "3.000", "h"),
            ("03:30:", "210.000",  "mn"),
            ("00:00:50", "50.000", "s")
        ]
        self.maccor_values_to_biologic_value_and_unit_test(
            converter._convert_time,
            tests,
        )

    def proc_step_to_seq_test(self, test_step_xml, diff_dict):
        """
        test utility for testing proc_step_to_seq
         """
        proc = xmltodict.parse(test_step_xml)
        test_step = proc["MaccorTestProcedure"]["ProcSteps"]["TestStep"]
        converter = MaccorToBiologicMb()

        expected = converter.blank_seq.copy()
        expected["Ns"] = 0
        expected["lim1_seq"] = 1
        expected["lim2_seq"] = 1
        expected["lim3_seq"] = 1
        expected.update(diff_dict)

        seq_num_by_step_num = {
            1: 0,
            2: 1,
        }

        result = converter._proc_step_to_seq(test_step, 1, seq_num_by_step_num, 0, 2)
        for key, value in expected.items():
            self.assertEqual(
                value,
                result[key],
                msg="Expected {0}: {1} got {0}: {2}".format(key, value, result[key]),
            )

    def test_rest_step_conversion(self):
        xml = (
            '<?xml version="1.0" encoding="UTF-8"?>'
            "<MaccorTestProcedure>"
            "  <ProcSteps>"
            "    <TestStep>"
            "      <StepType>  Rest  </StepType>"
            "      <StepMode>        </StepMode>"
            "      <StepValue></StepValue>"
            "      <Limits/>"
            "      <Ends>"
            "        <EndEntry>"
            "          <EndType>Voltage </EndType>"
            "          <SpecialType> </SpecialType>"
            "          <Oper>&gt;= </Oper>"
            "          <Step>002</Step>"
            "          <Value>4.4</Value>"
            "        </EndEntry>"
            "        <EndEntry>"
            "          <EndType>Voltage </EndType>"
            "          <SpecialType> </SpecialType>"
            "          <Oper>&lt;= </Oper>"
            "          <Step>002</Step>"
            "          <Value>2.5</Value>"
            "        </EndEntry>"
            "      </Ends>"
            "      <Reports>"
            "        <ReportEntry>"
            "          <ReportType>Voltage</ReportType>"
            "          <Value>2.2</Value>"
            "        </ReportEntry>"
            "      </Reports>"
            "      <Range>A</Range>"
            "      <Option1>N</Option1>"
            "      <Option2>N</Option2>"
            "      <Option3>N</Option3>"
            "      <StepNote></StepNote>"
            "    </TestStep>"
            "  </ProcSteps>"
            "</MaccorTestProcedure>"
        )
        diff_dict = {
            "ctrl_type": "Rest",
            "Apply I/C": "I",
            "N": "1.00",
            "charge/discharge": "Charge",
            "lim_nb": 2,
            "lim1_type": "Ecell",
            "lim1_comp": ">",
            "lim1_value": "4.400",
            "lim1_value_unit": "V",
            "lim2_type": "Ecell",
            "lim2_comp": "<",
            "lim2_value": "2.500",
            "lim2_value_unit": "V",
            "rec_nb": 1,
            "rec1_type": "Ecell",
            "rec1_value": "2.200",
            "rec1_value_unit": "V",
            "I Range": "10 A",
        }

        self.proc_step_to_seq_test(xml, diff_dict)
        pass

    def test_discharge_current_step_conversion(self):
        xml = (
            '<?xml version="1.0" encoding="UTF-8"?>'
            "<MaccorTestProcedure>"
            "  <ProcSteps>"
            "    <TestStep>"
            #      mispelling taken directly from sample file
            "      <StepType>Dischrge</StepType>"
            "      <StepMode>Current </StepMode>"
            "      <StepValue>1.0</StepValue>"
            "      <Limits/>"
            "      <Ends>"
            "        <EndEntry>"
            "          <EndType>StepTime</EndType>"
            "          <SpecialType> </SpecialType>"
            "          <Oper> = </Oper>"
            "          <Step>002</Step>"
            "          <Value>00:00:30</Value>"
            "        </EndEntry>"
            "        <EndEntry>"
            "          <EndType>Voltage </EndType>"
            "          <SpecialType> </SpecialType>"
            "          <Oper>&lt;= </Oper>"
            "          <Step>002</Step>"
            "          <Value>2.7</Value>"
            "        </EndEntry>"
            "        <EndEntry>"
            "          <EndType>Voltage </EndType>"
            "          <SpecialType> </SpecialType>"
            "          <Oper>&gt;= </Oper>"
            "          <Step>002</Step>"
            "          <Value>4.4</Value>"
            "        </EndEntry>"
            "      </Ends>"
            "      <Reports>"
            "        <ReportEntry>"
            "          <ReportType>Voltage </ReportType>"
            "          <Value>0.001</Value>"
            "        </ReportEntry>"
            "        <ReportEntry>"
            "          <ReportType>StepTime</ReportType>"
            #          10ms
            "          <Value>::.01</Value>"
            "        </ReportEntry>"
            "      </Reports>"
            "      <Range>A</Range>"
            "      <Option1>N</Option1>"
            "      <Option2>N</Option2>"
            "      <Option3>N</Option3>"
            "      <StepNote></StepNote>"
            "    </TestStep>"
            "  </ProcSteps>"
            "</MaccorTestProcedure>"
        )
        diff_dict = {
            "ctrl_type": "CC",
            "Apply I/C": "I",
            "ctrl1_val": "1.000",
            "ctrl1_val_unit": "A",
            "ctrl1_val_vs": "<None>",
            "N": "15.00",
            "charge/discharge": "Discharge",
            "lim_nb": 3,
            "lim1_type": "Time",
            "lim1_comp": ">",
            "lim1_value": "30.000",
            "lim1_value_unit": "s",
            "lim2_type": "Ecell",
            "lim2_comp": "<",
            "lim2_value": "2.700",
            "lim2_value_unit": "V",
            "lim3_type": "Ecell",
            "lim3_comp": ">",
            "lim3_value": "4.400",
            "lim3_value_unit": "V",
            "rec_nb": 2,
            "rec1_type": "Ecell",
            "rec1_value": "1.000",
            "rec1_value_unit": "mV",
            "rec2_type": "Time",
            "rec2_value": "10.000",
            "rec2_value_unit": "ms",
            "I Range": "10 A",
        }

        self.proc_step_to_seq_test(xml, diff_dict)
        pass

    def test_sample_conversion(self):
        # this test is super long but we want an E2E test that doesn't require us to inspect
        # a file and worry about things like encoding. Putting these strings in their own .py
        # file breaks our testing strategy so we're just keeping them here.

        maccor_ast = xmltodict.parse((
            "<MaccorTestProcedure>"
            "  <header>"
            "    <BuildTestVersion>"
            "      <major>1</major>"
            "      <minor>5</minor>"
            "      <release>7006</release>"
            "      <build>32043</build>"
            "    </BuildTestVersion>"
            "    <FileFormatVersion>"
            "      <BTVersion>11</BTVersion>"
            "    </FileFormatVersion>"
            "    <ProcDesc>"
            "      <desc></desc>"
            "    </ProcDesc>"
            "  </header>"
            "  <ProcSteps>"
            "    <TestStep>"
            "      <StepType>  Rest  </StepType>"
            "      <StepMode>        </StepMode>"
            "      <StepValue></StepValue>"
            "      <Limits/>"
            "      <Ends>"
            "        <EndEntry>"
            "          <EndType>StepTime</EndType>"
            "          <SpecialType> </SpecialType>"
            "          <Oper> = </Oper>"
            "          <Step>002</Step>"
            "          <Value>03:00:00</Value>"
            "        </EndEntry>"
            "      </Ends>"
            "      <Reports>"
            "        <ReportEntry>"
            "          <ReportType>StepTime</ReportType>"
            "          <Value>00:00:30</Value>"
            "        </ReportEntry>"
            "      </Reports>"
            "      <Range>4</Range>"
            "      <Option1>N</Option1>"
            "      <Option2>N</Option2>"
            "      <Option3>N</Option3>"
            "      <StepNote></StepNote>"
            "    </TestStep>"
            "    <TestStep>"
            "      <StepType> Charge </StepType>"
            "      <StepMode>Current </StepMode>"
            "      <StepValue>1.0</StepValue>"
            "      <Limits/>"
            "      <Ends>"
            "        <EndEntry>"
            "          <EndType>StepTime</EndType>"
            "          <SpecialType> </SpecialType>"
            "          <Oper> = </Oper>"
            "          <Step>003</Step>"
            "          <Value>00:00:01</Value>"
            "        </EndEntry>"
            "      </Ends>"
            "      <Reports>"
            "        <ReportEntry>"
            "          <ReportType>StepTime</ReportType>"
            "          <Value>::.01</Value>"
            "        </ReportEntry>"
            "      </Reports>"
            "      <Range>A</Range>"
            "      <Option1>N</Option1>"
            "      <Option2>N</Option2>"
            "      <Option3>N</Option3>"
            "      <StepNote>resistance check</StepNote>"
            "    </TestStep>"
            "    <TestStep>"
            "      <StepType>  Do 1  </StepType>"
            "      <StepMode>        </StepMode>"
            "      <StepValue></StepValue>"
            "      <Limits/>"
            "      <Ends/>"
            "      <Reports/>"
            "      <Range></Range>"
            "      <Option1></Option1>"
            "      <Option2></Option2>"
            "      <Option3></Option3>"
            "      <StepNote></StepNote>"
            "    </TestStep>"
            "    <TestStep>"
            "      <StepType> Charge </StepType>"
            "      <StepMode>Voltage </StepMode>"
            "      <StepValue>3.3</StepValue>"
            "      <Limits/>"
            "      <Ends>"
            "        <EndEntry>"
            "          <EndType>Current </EndType>"
            "          <SpecialType> </SpecialType>"
            "          <Oper>&lt;= </Oper>"
            "          <Step>005</Step>"
            "          <Value>0.0286</Value>"
            "        </EndEntry>"
            "      </Ends>"
            "      <Reports>"
            "        <ReportEntry>"
            "          <ReportType>Voltage </ReportType>"
            "          <Value>0.001</Value>"
            "        </ReportEntry>"
            "        <ReportEntry>"
            "          <ReportType>StepTime</ReportType>"
            "          <Value>00:02:00</Value>"
            "        </ReportEntry>"
            "      </Reports>"
            "      <Range>A</Range>"
            "      <Option1>N</Option1>"
            "      <Option2>N</Option2>"
            "      <Option3>N</Option3>"
            "      <StepNote>reset cycle C/20</StepNote>"
            "    </TestStep>"
            "    <TestStep>"
            "      <StepType> Loop 1 </StepType>"
            "      <StepMode>        </StepMode>"
            "      <StepValue></StepValue>"
            "      <Limits/>"
            "      <Ends>"
            "        <EndEntry>"
            "          <EndType>Loop Cnt</EndType>"
            "          <SpecialType> </SpecialType>"
            "          <Oper> = </Oper>"
            "          <Step>007</Step>"
            "          <Value>2</Value>"
            "        </EndEntry>"
            "      </Ends>"
            "      <Reports/>"
            "      <Range></Range>"
            "      <Option1></Option1>"
            "      <Option2></Option2>"
            "      <Option3></Option3>"
            "      <StepNote></StepNote>"
            "    </TestStep>"
            "    <TestStep>"
            "      <StepType>  Do 1  </StepType>"
            "      <StepMode>        </StepMode>"
            "      <StepValue></StepValue>"
            "      <Limits/>"
            "      <Ends/>"
            "      <Reports/>"
            "      <Range></Range>"
            "      <Option1></Option1>"
            "      <Option2></Option2>"
            "      <Option3></Option3>"
            "      <StepNote></StepNote>"
            "    </TestStep>"
            "    <TestStep>"
            "      <StepType> Charge </StepType>"
            "      <StepMode>Current </StepMode>"
            "      <StepValue>0.1429</StepValue>"
            "      <Limits>"
            "        <Voltage>4.2</Voltage>"
            "      </Limits>"
            "      <Ends>"
            "        <EndEntry>"
            "          <EndType>Current </EndType>"
            "          <SpecialType> </SpecialType>"
            "          <Oper>&lt;= </Oper>"
            "          <Step>008</Step>"
            "          <Value>0.0286</Value>"
            "        </EndEntry>"
            "      </Ends>"
            "      <Reports>"
            "        <ReportEntry>"
            "          <ReportType>Voltage </ReportType>"
            "          <Value>0.001</Value>"
            "        </ReportEntry>"
            "        <ReportEntry>"
            "          <ReportType>StepTime</ReportType>"
            "          <Value>00:02:00</Value>"
            "        </ReportEntry>"
            "      </Reports>"
            "      <Range>A</Range>"
            "      <Option1>N</Option1>"
            "      <Option2>N</Option2>"
            "      <Option3>N</Option3>"
            "      <StepNote>reset cycle C/20</StepNote>"
            "    </TestStep>"
            "    <TestStep>"
            "      <StepType>AdvCycle</StepType>"
            "      <StepMode>        </StepMode>"
            "      <StepValue></StepValue>"
            "      <Limits/>"
            "      <Ends/>"
            "      <Reports/>"
            "      <Range></Range>"
            "      <Option1></Option1>"
            "      <Option2></Option2>"
            "      <Option3></Option3>"
            "      <StepNote></StepNote>"
            "    </TestStep>"
            "    <TestStep>"
            "      <StepType> Loop 1 </StepType>"
            "      <StepMode>        </StepMode>"
            "      <StepValue></StepValue>"
            "      <Limits/>"
            "      <Ends>"
            "        <EndEntry>"
            "          <EndType>Loop Cnt</EndType>"
            "          <SpecialType> </SpecialType>"
            "          <Oper> = </Oper>"
            "          <Step>01</Step>"
            "          <Value>2</Value>"
            "        </EndEntry>"
            "      </Ends>"
            "      <Reports/>"
            "      <Range></Range>"
            "      <Option1></Option1>"
            "      <Option2></Option2>"
            "      <Option3></Option3>"
            "      <StepNote></StepNote>"
            "    </TestStep>"
            "    <TestStep>"
            "      <StepType>AdvCycle</StepType>"
            "      <StepMode>        </StepMode>"
            "      <StepValue></StepValue>"
            "      <Limits/>"
            "      <Ends/>"
            "      <Reports/>"
            "      <Range></Range>"
            "      <Option1></Option1>"
            "      <Option2></Option2>"
            "      <Option3></Option3>"
            "      <StepNote></StepNote>"
            "    </TestStep>"
            "    <TestStep>"
            "      <StepType>  End   </StepType>"
            "      <StepMode>        </StepMode>"
            "      <StepValue></StepValue>"
            "      <Limits/>"
            "      <Ends/>"
            "      <Reports/>"
            "      <Range></Range>"
            "      <Option1></Option1>"
            "      <Option2></Option2>"
            "      <Option3></Option3>"
            "      <StepNote></StepNote>"
            "    </TestStep>"
            "  </ProcSteps>"
            "</MaccorTestProcedure>"  
        ))

    

        expected_output = (
            "BT-LAB SETTING FILE\r\n"
            "\r\n"
            "Number of linked techniques : 1\r\n"
            "\r\n"
            "Filename : C:\\Users\\User\\Documents\\BT-Lab\\Data\\Grace\\BASF\\BCS - 171.64.160.115_Ja9_cOver70_CE3.mps\r\n\r\n"  # noqa
            "Device : BCS-805\r\n"
            "Ecell ctrl range : min = 0.00 V, max = 10.00 V\r\n"
            "Electrode material : \r\n"
            "Initial state : \r\n"
            "Electrolyte : \r\n"
            "Comments : \r\n"
            "Mass of active material : 0.001 mg\r\n"
            " at x = 0.000\r\n"  # leading space intentional
            "Molecular weight of active material (at x = 0) : 0.001 g/mol\r\n"
            "Atomic weight of intercalated ion : 0.001 g/mol\r\n"
            "Acquisition started at : xo = 0.000\r\n"
            "Number of e- transfered per intercalated ion : 1\r\n"
            "for DX = 1, DQ = 26.802 mA.h\r\n"
            "Battery capacity : 1.000 A.h\r\n"
            "Electrode surface area : 0.001 cm\N{superscript two}\r\n"
            "Characteristic mass : 8.624 mg\r\n"
            "Cycle Definition : Charge/Discharge alternance\r\n"
            "Do not turn to OCV between techniques\r\n"
            "\r\n"
            "Technique : 1\r\n"
            "Modulo Bat\r\n"
            "Ns                  0                   1                   2                   3                   4                   5                   6                   7                   8                   9                   \r\n"
            "ctrl_type           Rest                CC                  CV                  CV                  CV                  CC                  CV                  Loop                Loop                Loop                \r\n"
            "Apply I/C           I                   I                   I                   I                   I                   I                   I                   I                   I                   I                   \r\n"
            "ctrl1_val                               1.000               3.300               3.300               3.300               142.900             4.200               100.000             100.000             100.000             \r\n"
            "ctrl1_val_unit                          A                   V                   V                   V                   mA                  V                                                                               \r\n"
            "ctrl1_val_vs                            <None>              Ref                 Ref                 Ref                 <None>              Ref                                                                             \r\n"
            "ctrl2_val                                                                                                                                                                                                                   \r\n"
            "ctrl2_val_unit                                                                                                                                                                                                              \r\n"
            "ctrl2_val_vs                                                                                                                                                                                                                \r\n"
            "ctrl3_val                                                                                                                                                                                                                   \r\n"
            "ctrl3_val_unit                                                                                                                                                                                                              \r\n"
            "ctrl3_val_vs                                                                                                                                                                                                                \r\n"
            "N                   1.00                15.00               15.00               15.00               15.00               15.00               15.00                                                                           \r\n"
            "charge/discharge    Charge              Charge              Charge              Charge              Charge              Charge              Charge                                                                          \r\n"
            "ctrl_seq            0                   0                   0                   0                   0                   0                   0                   5                   0                   8                   \r\n"
            "ctrl_repeat         0                   0                   0                   0                   0                   0                   0                   1                   0                   1                   \r\n"
            "ctrl_trigger        Falling Edge        Falling Edge        Falling Edge        Falling Edge        Falling Edge        Falling Edge        Falling Edge        Falling Edge        Falling Edge        Falling Edge        \r\n"
            "ctrl_TO_t           0.000               0.000               0.000               0.000               0.000               0.000               0.000               0.000               0.000               0.000               \r\n"
            "ctrl_TO_t_unit      d                   d                   d                   d                   d                   d                   d                   d                   d                   d                   \r\n"
            "ctrl_Nd             6                   6                   6                   6                   6                   6                   6                   6                   6                   6                   \r\n"
            "ctrl_Na             1                   1                   1                   1                   1                   1                   1                   1                   1                   1                   \r\n"
            "ctrl_corr           1                   1                   1                   1                   1                   1                   1                   1                   1                   1                   \r\n"
            "lim_nb              1                   1                   1                   1                   1                   1                   1                   0                   0                   0                   \r\n"
            "lim1_type           Time                Time                |I|                 |I|                 |I|                 Ecell               |I|                 Time                Time                Time                \r\n"
            "lim1_comp           >                   >                   <                   <                   <                   >                   <                   <                   <                   <                   \r\n"
            "lim1_Q              Q limit             Q limit             Q limit             Q limit             Q limit             Q limit             Q limit             Q limit             Q limit             Q limit             \r\n"
            "lim1_value          3.000               1.000               28.600              28.600              28.600              4.200               28.600              0.000               0.000               0.000               \r\n"
            "lim1_value_unit     h                   s                   mA                  mA                  mA                  V                   mA                  s                   s                   s                   \r\n"
            "lim1_action         Next sequence       Next sequence       Next sequence       Next sequence       Next sequence       Next sequence       Next sequence       Next sequence       Next sequence       Next sequence       \r\n"
            "lim1_seq            1                   2                   3                   4                   5                   6                   7                   8                   9                   10                  \r\n"
            "lim2_type           Time                Time                Time                Time                Time                Time                Time                Time                Time                Time                \r\n"
            "lim2_comp           <                   <                   <                   <                   <                   <                   <                   <                   <                   <                   \r\n"
            "lim2_Q              Q limit             Q limit             Q limit             Q limit             Q limit             Q limit             Q limit             Q limit             Q limit             Q limit             \r\n"
            "lim2_value          0.000               0.000               0.000               0.000               0.000               0.000               0.000               0.000               0.000               0.000               \r\n"
            "lim2_value_unit     s                   s                   s                   s                   s                   s                   s                   s                   s                   s                   \r\n"
            "lim2_action         Next sequence       Next sequence       Next sequence       Next sequence       Next sequence       Next sequence       Next sequence       Next sequence       Next sequence       Next sequence       \r\n"
            "lim2_seq            1                   2                   3                   4                   5                   6                   7                   8                   9                   10                  \r\n"
            "lim3_type           Time                Time                Time                Time                Time                Time                Time                Time                Time                Time                \r\n"
            "lim3_comp           <                   <                   <                   <                   <                   <                   <                   <                   <                   <                   \r\n"
            "lim3_Q              Q limit             Q limit             Q limit             Q limit             Q limit             Q limit             Q limit             Q limit             Q limit             Q limit             \r\n"
            "lim3_value          0.000               0.000               0.000               0.000               0.000               0.000               0.000               0.000               0.000               0.000               \r\n"
            "lim3_value_unit     s                   s                   s                   s                   s                   s                   s                   s                   s                   s                   \r\n"
            "lim3_action         Next sequence       Next sequence       Next sequence       Next sequence       Next sequence       Next sequence       Next sequence       Next sequence       Next sequence       Next sequence       \r\n"
            "lim3_seq            1                   2                   3                   4                   5                   6                   7                   8                   9                   10                  \r\n"
            "rec_nb              1                   1                   2                   2                   2                   2                   2                   0                   0                   0                   \r\n"
            "rec1_type           Time                Time                Ecell               Ecell               Ecell               Ecell               Ecell               Time                Time                Time                \r\n"
            "rec1_value          30.000              10.000              1.000               1.000               1.000               1.000               1.000               10.000              10.000              10.000              \r\n"
            "rec1_value_unit     s                   ms                  mV                  mV                  mV                  mV                  mV                  s                   s                   s                   \r\n"
            "rec2_type           Time                Time                Time                Time                Time                Time                Time                Time                Time                Time                \r\n"
            "rec2_value          10.000              10.000              2.000               2.000               2.000               2.000               2.000               10.000              10.000              10.000              \r\n"
            "rec2_value_unit     s                   s                   mn                  mn                  mn                  mn                  mn                  s                   s                   s                   \r\n"
            "rec3_type           Time                Time                Time                Time                Time                Time                Time                Time                Time                Time                \r\n"
            "rec3_value          10.000              10.000              10.000              10.000              10.000              10.000              10.000              10.000              10.000              10.000              \r\n"
            "rec3_value_unit     s                   s                   s                   s                   s                   s                   s                   s                   s                   s                   \r\n"
            "E range min (V)     0.000               0.000               0.000               0.000               0.000               0.000               0.000               0.000               0.000               0.000               \r\n"
            "E range max (V)     10.000              10.000              10.000              10.000              10.000              10.000              10.000              10.000              10.000              10.000              \r\n"
            "I Range             10 A                10 A                10 A                10 A                10 A                10 A                10 A                1 mA                1 mA                1 mA                \r\n"
            "I Range min         Unset               Unset               Unset               Unset               Unset               Unset               Unset               Unset               Unset               Unset               \r\n"
            "I Range max         Unset               Unset               Unset               Unset               Unset               Unset               Unset               Unset               Unset               Unset               \r\n"
            "I Range init        Unset               Unset               Unset               Unset               Unset               Unset               Unset               Unset               Unset               Unset               \r\n"
            "auto rest           0                   0                   0                   0                   0                   0                   0                   0                   0                   0                   \r\n"
            "Bandwidth           4                   4                   4                   4                   4                   4                   4                   4                   4                   4                   \r\n"
        )

        expected_lines = expected_output.split("\r\n")

        converter = MaccorToBiologicMb()
        actual_output = converter.maccor_ast_to_protocol_str(maccor_ast, unroll=True, col_width=20)
        actual_lines = actual_output.split("\r\n")

        self.assertEqual(
           len(expected_lines),
           len(actual_lines),
        )    
        for i in range(0, len(expected_lines)):
            msg="At line {} expected:\n\"{}\"\ngot:\n\"{}\"".format(i + 1, expected_lines[i], actual_lines[i])

            self.assertEqual(
               expected_lines[i],
               actual_lines[i],
               msg
            )
        pass

    def test_remove_end_entries_by_pred(self):
        converter = MaccorToBiologicMb()
        fp = os.path.join(TEST_FILE_DIR, "diagnosticV4.000")
        maccor_ast = converter.load_maccor_ast(fp)

        def pred(end_entry, step_num):
            goto_step = int(end_entry["Step"])
            # filter all goto step 70s, except when that is Next Step
            return  goto_step != 70 or step_num == 69

        filtered_ast = converter.remove_end_entries_by_pred(maccor_ast, pred)
        filtered_steps = filtered_ast["MaccorTestProcedure"]["ProcSteps"]["TestStep"]
        
        step_1_filtered_end_entry = filtered_steps[0]["Ends"]["EndEntry"]
        self.assertEqual(
            OrderedDict,
            type(step_1_filtered_end_entry),
        )
        self.assertEqual(
            step_1_filtered_end_entry["Step"],
            "002"
        )

        step_68_filtered_end_entry = filtered_steps[68]["Ends"]["EndEntry"]
        self.assertEqual(
            OrderedDict,
            type(step_68_filtered_end_entry),
        )
        self.assertEqual(
            step_68_filtered_end_entry["Step"],
            "070"
        )
        pass

    def test_convert_diagnostic(self):
        # convert_diagnostic_v5_multi_techniques(source_file="BioTest_000001.000")
        with ScratchDir("."):
            convert_diagnostic_v5_multi_techniques(source_file="diagnosticV5.000")
