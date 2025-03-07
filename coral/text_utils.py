import re

COMMA_NUMBER_REGEX = re.compile(r"\d{1,3}(?:,\d{3})*")
AMPLICON_SEPARATOR = "-----------------------------------------------"
SUMMARY_HEADER = "CoRAL v"
