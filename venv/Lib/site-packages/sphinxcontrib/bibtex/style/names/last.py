from pybtex.style.names import BaseNameStyle, name_part
from pybtex.style.template import join


class LastNameStyle(BaseNameStyle):
    """A simple name style for formatting the last name of an author."""

    def format(self, person, abbr=True):
        """Format last names."""

        return join[
            name_part(tie=True)[person.rich_prelast_names],
            name_part[person.rich_last_names],
        ]
