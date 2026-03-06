from pydantic import BaseSettings, Field, validator


class ExampleSettings(BaseSettings):
    """Document your project settings very conveniently. Applies like wise
    to pydantic models.

    """

    field_plain_with_validator: int = 100
    """Show standard field with type annotation."""

    field_plain = 100
    """Show standard field without type annotation."""

    field_with_validator_and_alias: str = Field("FooBar", alias="BarFoo", env="BarFoo")
    """Shows corresponding validator with link/anchor."""

    field_with_validator_and_alias_typeless = Field(
        "FooBar", alias="BarFoo", env="BarFoo"
    )
    """Shows corresponding validator with link/anchor but without type annotation."""

    field_with_constraints_and_description: int = Field(
        default=5, ge=0, le=100, description="Shows constraints within doc string."
    )

    @validator("field_with_validator_and_alias", "field_plain_with_validator")
    def check_max_length_ten(cls, v):
        """Show corresponding field with link/anchor."""

        if not len(v) < 10:
            raise ValueError("No more than 10 characters allowed")

    class Config:
        env_prefix = "foo_"
        allow_mutation = True

    def method_without_sig(self) -> BaseSettings:
        return BaseSettings()
